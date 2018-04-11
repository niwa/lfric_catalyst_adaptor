#include <string>
#include <mpi.h>

#include <vtkCPAdaptorAPI.h>
#include <vtkCPProcessor.h>
#include <vtkCPPipeline.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkSmartPointer.h>

#include "vtk_pipeline.h"

//
// Define wrappers for coprocessor interface
//

// We use Catalyst's "CPAdaptorAPI" layer here, which is a C wrapper
// API for talking to Catalyst, register our visualisation pipeline etc.
// Users can choose between a C++ pipeline and a scripted Python pipeline.
//
// The adaptor may need to be reimplemented with the C++ API for more
// fine-grained control over Catalyst.

extern "C" {

  void coprocessor_initialize(const int visualisationFrequency, const char * outputFileName,
                              const MPI_Fint mpi_fortran_comm, const int usePythonPipeline,
                              const char * pythonScript) {

    // Check if the coprocessor has already been initialised, only clear out pipelines in that case
    if (!vtkCPAdaptorAPI::GetCoProcessor()) {
      // Note that this will initialise Catalyst with MPI_COMM_WORLD communicator. This seems ok
      // at the moment, but may need to be revisited. The Catalyst C wrapper API does not expose
      // a method to initialise Catalyst with a different communicator, but the C++ API has an
      // "initialize" method for this purpose.
      vtkCPAdaptorAPI::CoProcessorInitialize();
    }
    else {
      vtkCPAdaptorAPI::GetCoProcessor()->RemoveAllPipelines();
    }

    if (usePythonPipeline) {

      if (!strcmp(pythonScript,"")) {
	vtkGenericWarningMacro("coprocessor_initialize: No filename provided for pythonScript.");
        return;
      }

      // Create new Python pipeline object
      vtkSmartPointer<vtkCPPythonScriptPipeline> pipeline = vtkSmartPointer<vtkCPPythonScriptPipeline>::New();
      pipeline->Initialize(pythonScript);
      vtkCPAdaptorAPI::GetCoProcessor()->AddPipeline(pipeline);

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

      // Retrieve C communicator from Fortran
      MPI_Comm mpi_comm = MPI_Comm_f2c(mpi_fortran_comm);

      // Find out MPI rank and size to handle grid partitioning in VTK grid
      int mpiSize = 1;
      int mpiRank = 0;
      MPI_Comm_rank(mpi_comm, &mpiRank);
      MPI_Comm_size(mpi_comm, &mpiSize);

      // Create new C++ pipeline object
      vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();
      std::string outputFileNameString(outputFileName);
      pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
      vtkCPAdaptorAPI::GetCoProcessor()->AddPipeline(pipeline);

    }

  }

  void coprocessor_requestdatadescription(int * timeStep, double * time, int * coprocessThisTimeStep) {
    // Note that a mechanism for forcing output can be implemented here using
    // vtkCPAdaptorAPI::GetCoProcessorData()->SetForceOutput(true);
    // before calling "RequestDataDescription".
    vtkCPAdaptorAPI::RequestDataDescription(timeStep, time, coprocessThisTimeStep);
  }

  void coprocessor_needtocreategrid(int * needGrid) {
    vtkCPAdaptorAPI::NeedToCreateGrid(needGrid);
  }

  void coprocessor_coprocess() {
    vtkCPAdaptorAPI::CoProcess();
  }

  void coprocessor_finalize() {
    vtkCPAdaptorAPI::CoProcessorFinalize();
  }

}
