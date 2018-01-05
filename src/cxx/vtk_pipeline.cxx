#include <string>
#include <sstream>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkUnstructuredGrid.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <vtkPVTrivialProducer.h>
#include <vtkCompleteArrays.h>
#include <vtkXMLPUnstructuredGridWriter.h>

#include "vtk_pipeline.h"

//
// Create visualisation pipeline
//

// The pipeline is implemented as a VTK object and includes a number of
// callback functions that will be used by the Catalyst coprocessor.

// Insert standard implementation of "New"
vtkStandardNewMacro(vtkCPVTKPipeline);

// Constructor and desctructor
vtkCPVTKPipeline::vtkCPVTKPipeline(){this->OutputFrequency = 0;}
vtkCPVTKPipeline::~vtkCPVTKPipeline(){}

// Set basic pipeline parameters
void vtkCPVTKPipeline::SetVTKPipelineParameters(int outputFrequency, std::string& fileName) {
  this->OutputFrequency = outputFrequency;
  this->FileName = fileName;
}

// Callback function: work out if we should produce output for this call
// and let Catalyst know
int vtkCPVTKPipeline::RequestDataDescription(vtkCPDataDescription* dataDescription) {

  if(!dataDescription) {
    vtkWarningMacro("vtk_pipeline: RequestDataDescription: dataDescription object is undefined.");
    return 0;
  }

  if(this->FileName.empty()) {
    vtkWarningMacro("vtk_pipeline: RequestDataDescription: No output filename was set.");
    return 0;
  }

  if(dataDescription->GetNumberOfInputDescriptions() != 1) {
    vtkWarningMacro("vtk_pipeline: RequestDataDescription: Expected exactly 1 input description.");
    return 0;
  }

  // Check if we need to produce output (forced or via request)
  if(dataDescription->GetForceOutput() == true ||
    (this->OutputFrequency != 0 && dataDescription->GetTimeStep() % this->OutputFrequency == 0) ) {
    dataDescription->GetInputDescription(0)->AllFieldsOn();
    dataDescription->GetInputDescription(0)->GenerateMeshOn();
    return 1;
  }
  return 0;
}

// Callback function
int vtkCPVTKPipeline::CoProcess(vtkCPDataDescription* dataDescription) {

  if(!dataDescription) {
    vtkWarningMacro("vtk_pipeline: CoProcess: DataDescription is NULL");
    return 0;
  }

  // Try to get grid
  vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(
    dataDescription->GetInputDescriptionByName("input")->GetGrid());
  if(!grid) {
    vtkWarningMacro("vtk_pipeline: CoProcess: DataDescription is missing input unstructured grid.");
    return 0;
  }

  vtkNew<vtkPVTrivialProducer> producer;
  producer->SetOutput(grid);

  // If process 0 doesn't have any points or cells, the writer may
  // have problems in parallel so we use completeArrays to fill in
  // the missing information.
  vtkNew<vtkCompleteArrays> completeArrays;
  completeArrays->SetInputConnection(producer->GetOutputPort());

  vtkNew<vtkXMLPUnstructuredGridWriter> writer;
  writer->SetInputConnection(completeArrays->GetOutputPort());
  std::ostringstream o;
  o << dataDescription->GetTimeStep();
  std::string name = this->FileName + o.str() + ".pvtu";
  writer->SetFileName(name.c_str());
  writer->Update();

  return 1;
}

void vtkCPVTKPipeline::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "OutputFrequency: " << this->OutputFrequency << "\n";
  os << indent << "FileName: " << this->FileName << "\n";
}
