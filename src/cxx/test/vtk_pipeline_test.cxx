#include "Catch2/catch.hpp"
#include "vtk_pipeline.h"

#include "vtkSmartPointer.h"
#include "vtkCPDataDescription.h"

TEST_CASE( "VTK Pipeline New method works", "[vtk_pipeline]" ) {

  vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();

  REQUIRE( pipeline != nullptr);
  REQUIRE( pipeline->GetVTKPipelineOutputFrequency() == 0 );
  REQUIRE( pipeline->GetVTKPipelineFileName().compare("catalyst_output.vtk") == 0 );
  REQUIRE( pipeline->GetVTKPipelineMPIRank() == 0 );
  REQUIRE( pipeline->GetVTKPipelineMPISize() == 1 );

}

TEST_CASE( "VTK Pipeline SetVTKPipelineParameters method works", "[vtk_pipeline]" ) {

  vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();

  const int visualisationFrequency = 1;
  const std::string outputFileNameString("testoutput.vtk");
  const int mpiRank = 2;
  const int mpiSize = 3;

  pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);

  REQUIRE( pipeline->GetVTKPipelineOutputFrequency() == 1 );
  REQUIRE( pipeline->GetVTKPipelineFileName().compare("testoutput.vtk") == 0 );
  REQUIRE( pipeline->GetVTKPipelineMPIRank() == 2 );
  REQUIRE( pipeline->GetVTKPipelineMPISize() == 3 );

}

TEST_CASE( "VTK Pipeline RequestDataDescription method works", "[vtk_pipeline]" ) {

  vtkSmartPointer<vtkCPDataDescription> testDataDescription = vtkSmartPointer<vtkCPDataDescription>::New();
  testDataDescription->AddInput("input");

  vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();

  int visualisationFrequency = 1;
  std::string outputFileNameString("output.vtk");
  const int mpiRank = 2;
  const int mpiSize = 3;

  SECTION( "No dataDescription provided" ) {
    REQUIRE( pipeline->RequestDataDescription(nullptr) == 0 );
  }

  SECTION( "No filename provided" ) {
    outputFileNameString.clear();
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 0 );
  }

  SECTION( "Too many inputs available" ) {
    testDataDescription->AddInput("another_input");
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 0 );
  }

  SECTION( "Force Output is on with nonsensical timeStep" ) {
    testDataDescription->SetForceOutput(true);
    double time = -1.0;
    int timeStep = -1;
    testDataDescription->SetTimeData(time, timeStep);
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 1 );
  }

  SECTION( "Force Output is off with nonsensical timeStep" ) {
    testDataDescription->SetForceOutput(false);
    double time = -1.0;
    int timeStep = -1;
    testDataDescription->SetTimeData(time, timeStep);
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 0 );
  }

  SECTION( "No output with visualisationFrequency 0" ) {
    testDataDescription->SetForceOutput(false);
    double time = 0.0;
    int timeStep = 0;
    testDataDescription->SetTimeData(time, timeStep);
    visualisationFrequency = 0;
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 0 );
  }

  SECTION( "Output at timeStep 0 with visualisationFrequency 1" ) {
    testDataDescription->SetForceOutput(false);
    double time = 0.0;
    int timeStep = 0;
    testDataDescription->SetTimeData(time, timeStep);
    visualisationFrequency = 1;
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 1 );
  }

  SECTION( "Output at timeStep 7 with visualisationFrequency 7" ) {
    testDataDescription->SetForceOutput(false);
    double time = 0.0;
    int timeStep = 7;
    testDataDescription->SetTimeData(time, timeStep);
    visualisationFrequency = 7;
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 1 );
  }

  SECTION( "No output at timeStep 8 with visualisationFrequency 7" ) {
    testDataDescription->SetForceOutput(false);
    double time = 0.0;
    int timeStep = 8;
    testDataDescription->SetTimeData(time, timeStep);
    visualisationFrequency = 7;
    pipeline->SetVTKPipelineParameters(visualisationFrequency, outputFileNameString, mpiRank, mpiSize);
    REQUIRE( pipeline->RequestDataDescription(testDataDescription) == 0 );
  }

}

TEST_CASE( "VTK Pipeline CoProcess method fail tests", "[vtk_pipeline]" ) {

  vtkSmartPointer<vtkCPVTKPipeline> pipeline = vtkSmartPointer<vtkCPVTKPipeline>::New();

  SECTION( "No dataDescription provided" ) {
    REQUIRE( pipeline->CoProcess(nullptr) == 0 );
  }

  SECTION( "No InputDescription named input set" ) {
    vtkSmartPointer<vtkCPDataDescription> testDataDescription = vtkSmartPointer<vtkCPDataDescription>::New();
    testDataDescription->AddInput("another_input");
    REQUIRE( pipeline->CoProcess(testDataDescription) == 0 );
  }

  SECTION( "No grid set" ) {
    vtkSmartPointer<vtkCPDataDescription> testDataDescription = vtkSmartPointer<vtkCPDataDescription>::New();
    testDataDescription->AddInput("input");
    REQUIRE( pipeline->CoProcess(testDataDescription) == 0 );
  }

}
