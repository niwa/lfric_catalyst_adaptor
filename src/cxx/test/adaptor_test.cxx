#include "Catch2/catch.hpp"
#include "adaptor.h"
#include "coprocessor.h"

#include <vtkCPInputDataDescription.h>
#include <vtkCellType.h>

// Class for defining test grids using points and connectivities
class TestGridDef {
  public:
    long npoints;
    long ncells;
    double * coords;
    long * cell_points;
    short * ghost_mask;

    TestGridDef (long new_npoints, long new_ncells, double new_coords[],
                 long new_cell_points[], short new_ghost_mask[]) {
      npoints = new_npoints;
      ncells = new_ncells;

      coords = new double[npoints*3];
      for (long i = 0; i < npoints*3; i++) {coords[i] = new_coords[i];}

      cell_points = new long[npoints];
      for (long i = 0; i < npoints; i++) {cell_points[i] = new_cell_points[i];}

      ghost_mask = new short[ncells];
      for (long i = 0; i < ncells; i++) {ghost_mask[i] = new_ghost_mask[i];}
    }

    ~TestGridDef () {
      delete coords;
      delete cell_points;
      delete ghost_mask;
    }
};

// Return a simple grid definition object with one cell
TestGridDef * NewSingleCellGrid() {
  double coords[] = {0.0, 0.0, 0.0,
                     1.0, 0.0, 0.0,
                     1.0, 1.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0,
                     1.0, 0.0, 1.0,
                     1.0, 1.0, 1.0,
                     0.0, 1.0, 1.0};
  long cell_points[] = {0, 1, 2, 3, 4, 5, 6, 7};
  short ghost_mask[] = {0};
  TestGridDef * griddef = new TestGridDef(8, 1, coords, cell_points, ghost_mask);
  return griddef;
}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Adaptor CreateGrid with single cell grid", "[coprocessor]" ) {

  // Keep the original dataDescription object and
  // create a second one for this test
  vtkCPDataDescription * dataDescription_save = dataDescription;
  dataDescription = vtkCPDataDescription::New();
  dataDescription->AddInput("input");

  // Get a simple grid with a single cell
  TestGridDef * griddef = NewSingleCellGrid();
  short use_ghost_mask = 0;
  short mirror_periodic = 0;

  // Call CreateGrid function
  adaptor_creategrid(griddef->coords, griddef->npoints, griddef->cell_points, griddef->ncells,
                     griddef->ghost_mask, use_ghost_mask, mirror_periodic);

  delete griddef;

  vtkCPInputDataDescription * InputDescription = dataDescription->GetInputDescription(0);
  vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(InputDescription->GetGrid());
  REQUIRE( grid != NULL );
  REQUIRE( grid->GetNumberOfCells() == 1 );
  REQUIRE( grid->GetCellType(0) == VTK_HEXAHEDRON );
  REQUIRE( grid->HasAnyGhostCells() == 0 );

  vtkPoints * gridPoints = grid->GetPoints();
  REQUIRE( gridPoints->GetNumberOfPoints() == 8 );

  double gridBounds[6];
  grid->GetBounds(gridBounds);
  REQUIRE( gridBounds[0] == Approx(0.0) );
  REQUIRE( gridBounds[1] == Approx(1.0) );
  REQUIRE( gridBounds[2] == Approx(0.0) );
  REQUIRE( gridBounds[3] == Approx(1.0) );
  REQUIRE( gridBounds[4] == Approx(0.0) );
  REQUIRE( gridBounds[5] == Approx(1.0) );

  // Clean up test object and restore original one
  dataDescription->Delete();
  dataDescription = dataDescription_save;
  dataDescription_save = NULL;

}
