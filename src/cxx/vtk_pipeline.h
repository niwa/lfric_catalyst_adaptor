#ifndef VTK_PIPELINE_H
#define VTK_PIPELINE_H

#include <vtkCPPipeline.h>
#include <string>

// This class extends the generic vtkCPPipeline interface
// class to define the coprocessing pipeline, create the
// logic that decides if the coprocessor should run, and
// additional coprocessor configuration.
//
// Class definition follows the C++ pipeline example in the
// Catalyst source code
class vtkCPVTKPipeline : public vtkCPPipeline {
 public:

  // Standard VTK members
  static vtkCPVTKPipeline* New();
  virtual void PrintSelf(ostream& os, vtkIndent indent);
  vtkTypeMacro(vtkCPVTKPipeline,vtkCPPipeline);

  // Our class implements these two Catalyst coprocessor API functions
  virtual int RequestDataDescription(vtkCPDataDescription* dataDescription);
  virtual int CoProcess(vtkCPDataDescription* dataDescription);

  // Set basic parameters of the visualisation pipeline
  virtual void SetVTKPipelineParameters(int outputFrequency, std::string& fileName);

 protected:
  vtkCPVTKPipeline();
  virtual ~vtkCPVTKPipeline();

 private:

  // Disable these (not implemented)
  vtkCPVTKPipeline(const vtkCPVTKPipeline&);
  void operator=(const vtkCPVTKPipeline&);

  int OutputFrequency;
  std::string FileName;
};
#endif
