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

// Catalyst provides generic API functions for this, but we need to set
// visualisation parameters and register our pipeline, so it seems
// best to do everything directly in C++ and provide C-style wrappers.
// It also enables users to choose between a C++ pipeline and a scripted
// Python pipeline.

extern "C" {

  void coprocessor_initialize(const int visualisationFrequency, const char * outputFileName,
                              const MPI_Comm mpi_comm, const int usePythonPipeline,
                              const char * pythonScript) {

    // Check if the coprocessor has already been initialised, only clear out pipelines in that case
    if (!vtkCPAdaptorAPI::GetCoProcessor()) {
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
    vtkCPAdaptorAPI::RequestDataDescription(timeStep, time, coprocessThisTimeStep);
  }

  void coprocessor_needtocreategrid(int * needGrid) {
    vtkCPAdaptorAPI::NeedToCreateGrid(needGrid);
  }

  void coprocessor_coprocess() {
    vtkCPAdaptorAPI::CoProcess();
  }

  void coprocessor_finalize() {
    vtkCPAdaptorAPI::CoProcessorFinalize ();
  }

}
