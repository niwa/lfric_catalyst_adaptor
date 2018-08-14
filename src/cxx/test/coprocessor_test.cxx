#include "Catch2/catch.hpp"
#include "coprocessor.h"

#include <vtkCPInputDataDescription.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

TEST_CASE( "C++ coprocessor initialised", "[coprocessor]" ) {

  // Objects must be allocated if initialisation was successful
  REQUIRE ( Processor != nullptr );
  REQUIRE ( dataDescription != nullptr );

}

TEST_CASE( "C++ coprocessor RequestDataDescription works", "[coprocessor]" ) {

  int coprocessThisTimeStep, timeStep;
  double time = 0.0;

  SECTION( "Output at time step -1 NOT required" ) {
    coprocessThisTimeStep = -1;
    timeStep = -1;
    coprocessor_requestdatadescription(&timeStep, &time, &coprocessThisTimeStep);
    REQUIRE ( coprocessThisTimeStep == 0 );
  }

  SECTION( "Output at time step 0 required" ) {
    coprocessThisTimeStep = -1;
    timeStep = 0;
    coprocessor_requestdatadescription(&timeStep, &time, &coprocessThisTimeStep);
    REQUIRE ( coprocessThisTimeStep == 1 );
  }

  SECTION( "Output at time step 1 NOT required" ) {
    coprocessThisTimeStep = -1;
    timeStep = 1;
    coprocessor_requestdatadescription(&timeStep, &time, &coprocessThisTimeStep);
    REQUIRE ( coprocessThisTimeStep == 0 );
  }

  SECTION( "Output at time step 2 required" ) {
    coprocessThisTimeStep = -1;
    timeStep = 2;
    coprocessor_requestdatadescription(&timeStep, &time, &coprocessThisTimeStep);
    REQUIRE ( coprocessThisTimeStep == 1 );
  }

}

TEST_CASE( "C++ coprocessor NeedToCreateGrid works", "[coprocessor]" ) {

  // Keep the original dataDescription object and
  // create a second one for this test
  vtkCPDataDescription * dataDescription_save = dataDescription;
  dataDescription = vtkCPDataDescription::New();
  dataDescription->AddInput("input");

  int needGrid;

  SECTION( "Grid required" ) {
    needGrid = -1;
    coprocessor_needtocreategrid(&needGrid);
    REQUIRE( needGrid == 1 );
  }

  SECTION( "Grid NOT required" ) {
    // Create dummy grid for this test
    vtkSmartPointer<vtkUnstructuredGrid> dummygrid = 
       vtkSmartPointer<vtkUnstructuredGrid>::New();
    dataDescription->GetInputDescription(0)->SetGrid(dummygrid);
    needGrid = -1;
    coprocessor_needtocreategrid(&needGrid);
    REQUIRE( needGrid == 0 );
  }

  // Clean up test object and restore original one
  dataDescription->Delete();
  dataDescription = dataDescription_save;
  dataDescription_save = nullptr;

}
