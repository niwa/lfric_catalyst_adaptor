#include "coprocessor.h"
#include "vtk_pipeline.h"

#include <vtkMPI.h>
#include <vtkCPPipeline.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkCPInputDataDescription.h>
#include <vtkSmartPointer.h>

#include <string>

vtkCPProcessor * Processor = NULL;
vtkCPDataDescription * dataDescription = NULL;

extern "C" {

  void coprocessor_initialize(const int visualisationFrequency, const char * outputFileName,
                              const MPI_Fint mpi_fortran_comm, const int usePythonPipeline,
                              const char * pythonScript) {

    // Retrieve C communicator from Fortran
    MPI_Comm mpi_comm = MPI_Comm_f2c(mpi_fortran_comm);

    // Check if the coprocessor has already been initialised, only clear out pipelines in that case
    if(!Processor) {
      Processor = vtkCPProcessor::New();
      // Initialise with custom MPI communicator, Catalyst assumes MPI_COMM_WORLD otherwise,
      // which will cause problems if mpi_comm is a split communicator
      vtkMPICommunicatorOpaqueComm mpi_vtk_comm = vtkMPICommunicatorOpaqueComm(&mpi_comm);
      Processor->Initialize(mpi_vtk_comm);
    }
    else {
      Processor->RemoveAllPipelines();
    }

    if (usePythonPipeline) {

      if (!strcmp(pythonScript,"")) {
	vtkGenericWarningMacro("coprocessor_initialize: No filename provided for pythonScript.");
        return;
      }

      // Create new Python pipeline object
      vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
      pipeline->Initialize(pythonScript);
      Processor->AddPipeline(pipeline);

    }
    else {

      if (visualisationFrequency < 0) {
	vtkGenericWarningMacro("coprocessor_initialize: Invalid number provided for visualisationFrequency:" << visualisationFrequency);
	return;
      }
      if (!strcmp(outputFileName,"")) {
	vtkGenericWarningMacro("coprocessor_initialize: No filename provided for outputFileName.");
	return;
      }

      // Find out MPI rank and size to handle grid partitioning in VTK grid
      int mpiSize = 1;
      int mpiRank = 0;
      MPI_Comm_rank(mpi_comm, &mpiRank);
      MPI_Comm_size(mpi_comm, &mpiSize);

      // Create new C++ pipeline object
      vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();
      std::string outputFileNameString(outputFileName);
      pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
      Processor->AddPipeline(pipeline);

    }

    // Create DataDescription object that keeps track of simulation data
    // Also add new InputDataDescription for our VTK grid and data structures,
    // using default name "input"
    if (!dataDescription) {
      dataDescription = vtkCPDataDescription::New();
      dataDescription->AddInput("input");
    }

  }

  void coprocessor_requestdatadescription(int * timeStep, double * time, int * coprocessThisTimeStep) {
    if (!Processor) {
      vtkGenericWarningMacro("coprocessor_requestdatadescription: Unable to access Processor.");
      return;
    }
    if (!dataDescription) {
      vtkGenericWarningMacro("coprocessor_requestdatadescription: Unable to access dataDescription.");
      return;
    }
    // Check if anything needs to be done for this simulation time and timeStep
    dataDescription->SetTimeData(*time, *timeStep);
    *coprocessThisTimeStep = Processor->RequestDataDescription(dataDescription);
  }

  void coprocessor_needtocreategrid(int * needGrid) {
    if (!Processor) {
      vtkGenericWarningMacro("coprocessor_needtocreategrid: Unable to access Processor.");
      return;
    }
    if (!dataDescription) {
      vtkGenericWarningMacro("coprocessor_needtocreategrid: Unable to access dataDescription.");
      return;
    }
    if (dataDescription->GetNumberOfInputDescriptions() != 1) {
      vtkGenericWarningMacro("coprocessor_needtocreategrid: expected exactly 1 input description.");
      return;
    }
    // Check if a grid has already been registered with our inputDescription
    if (!dataDescription->GetInputDescription(0)->GetGrid()) {
      *needGrid = 1;
    } else {
      *needGrid = 0;
    }
  }

  void coprocessor_coprocess() {
    if (!Processor) {
      vtkGenericWarningMacro("coprocessor_coprocess: Unable to access Processor.");
      return;
    }
    if (!dataDescription) {
      vtkGenericWarningMacro("coprocessor_coprocess: Unable to access dataDescription.");
      return;
    }
    int retval = Processor->CoProcess(dataDescription);
    if (retval != 1) {
      vtkGenericWarningMacro("coprocessor_coprocess: Processor reported failure.");
    }
  }

  void coprocessor_finalize() {
    if (Processor) {
      Processor->Finalize();
      Processor->Delete();
      Processor = NULL;
    }
    if (dataDescription) {
      dataDescription->Delete();
      dataDescription = NULL;
    }
  }

}
