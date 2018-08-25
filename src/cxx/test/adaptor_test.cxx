#include "Catch2/catch.hpp"
#include "adaptor.h"
#include "coprocessor.h"

#include <vtkCPInputDataDescription.h>
#include <vtkCellType.h>
#include <vtkDoubleArray.h>

// Class for defining test grids using points and connectivities
class TestGridDef {
  public:
    long npoints;
    long ncells;
    double * coords;
    long * cell_points;
    short * ghost_mask;
    const long npoints_per_cell = 8; // Hexahedron

    TestGridDef (long new_npoints, long new_ncells, double new_coords[],
                 long new_cell_points[], short new_ghost_mask[]) {
      npoints = new_npoints;
      ncells = new_ncells;

      coords = new double[npoints*3];
      for (long i = 0; i < npoints*3; i++) {coords[i] = new_coords[i];}

      cell_points = new long[ncells*npoints_per_cell];
      for (long i = 0; i < ncells*npoints_per_cell; i++) {cell_points[i] = new_cell_points[i];}

      ghost_mask = new short[ncells];
      for (long i = 0; i < ncells; i++) {ghost_mask[i] = new_ghost_mask[i];}
    }

    ~TestGridDef () {
      delete[] coords;
      delete[] cell_points;
      delete[] ghost_mask;
    }

};

// Return a simple grid definition object with one cell
// marked as a ghost cell
TestGridDef NewSingleCellGrid() {
  double coords[] = {0.0, 0.0, 0.0,
                     1.0, 0.0, 0.0,
                     1.0, 1.0, 0.0,
                     0.0, 1.0, 0.0,
                     0.0, 0.0, 1.0,
                     1.0, 0.0, 1.0,
                     1.0, 1.0, 1.0,
                     0.0, 1.0, 1.0};
  long cell_points[] = {0, 1, 2, 3, 4, 5, 6, 7};
  short ghost_mask[] = {1};
  return TestGridDef(8, 1, coords, cell_points, ghost_mask);
}

// Six cells in the xy plane, periodic in x direction
TestGridDef NewSixCellPeriodicGrid() {
  double coords[] = {0.0, 0.0, 0.0,
                     1.0, 0.0, 0.0,
                     2.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     1.0, 1.0, 0.0,
                     2.0, 1.0, 0.0,
                     0.0, 2.0, 0.0,
                     1.0, 2.0, 0.0,
                     2.0, 2.0, 0.0,
                     0.0, 0.0, 1.0,
                     1.0, 0.0, 1.0,
                     2.0, 0.0, 1.0,
                     0.0, 1.0, 1.0,
                     1.0, 1.0, 1.0,
                     2.0, 1.0, 1.0,
                     0.0, 2.0, 1.0,
                     1.0, 2.0, 1.0,
                     2.0, 2.0, 1.0};
  long cell_points[] = {0, 1, 4, 3, 9, 10, 13, 12, // South West
                        1, 2, 5, 4, 10, 11, 14, 13, // South East
                        3, 4, 7, 6, 12, 13, 16, 15, // North West
                        4, 5, 8, 7, 13, 14, 17, 16, // North East
                        2, 0, 3, 5, 11, 9, 12, 14, // Periodic SE
                        5, 3, 6, 8, 14, 12, 15, 17}; // Periodic NE
  short ghost_mask[] = {0, 0, 0, 0, 0, 0};
  return TestGridDef(18, 6, coords, cell_points, ghost_mask);
}

// Nine cells in the xy plane, periodic in x and y directions
TestGridDef NewNineCellBiPeriodicGrid() {
  double coords[] = {0.0, 0.0, 0.0,
                     1.0, 0.0, 0.0,
                     2.0, 0.0, 0.0,
                     0.0, 1.0, 0.0,
                     1.0, 1.0, 0.0,
                     2.0, 1.0, 0.0,
                     0.0, 2.0, 0.0,
                     1.0, 2.0, 0.0,
                     2.0, 2.0, 0.0,
                     0.0, 0.0, 1.0,
                     1.0, 0.0, 1.0,
                     2.0, 0.0, 1.0,
                     0.0, 1.0, 1.0,
                     1.0, 1.0, 1.0,
                     2.0, 1.0, 1.0,
                     0.0, 2.0, 1.0,
                     1.0, 2.0, 1.0,
                     2.0, 2.0, 1.0};
  long cell_points[] = {0, 1, 4, 3, 9, 10, 13, 12, // South West
                        1, 2, 5, 4, 10, 11, 14, 13, // South East
                        3, 4, 7, 6, 12, 13, 16, 15, // North West
                        4, 5, 8, 7, 13, 14, 17, 16, // North East
                        2, 0, 3, 5, 11, 9, 12, 14, // Periodic SE
                        5, 3, 6, 8, 14, 12, 15, 17, // Periodic NE
                        8, 6, 0, 2, 17, 15, 9, 11, // Periodic corner
                        6, 7, 1, 0, 15, 16, 10, 9, // Periodic NW
                        7, 8, 2, 1, 16, 17, 11, 10}; // Periodic NE
  short ghost_mask[] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
  return TestGridDef(18, 9, coords, cell_points, ghost_mask);
}

// ------------------------------------------------------------------------------------------

TEST_CASE( "Adaptor CreateGrid with single cell grid", "[coprocessor]" ) {

  // Keep the original dataDescription object and
  // create a second one for this test
  vtkCPDataDescription * dataDescription_save = dataDescription;
  dataDescription = vtkCPDataDescription::New();
  dataDescription->AddInput("input");

  // Get a simple grid with a single cell
  TestGridDef griddef = NewSingleCellGrid();
  short mirror_periodic = 0;

  SECTION( "Basic tests") {

    short use_ghost_mask = 0;

    // Call CreateGrid function
    adaptor_creategrid(griddef.coords, griddef.npoints, griddef.cell_points, griddef.ncells,
                       griddef.ghost_mask, use_ghost_mask, mirror_periodic);

    vtkCPInputDataDescription * InputDescription = dataDescription->GetInputDescription(0);
    vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(InputDescription->GetGrid());
    REQUIRE( grid != nullptr );
    REQUIRE( grid->IsTypeOf("vtkUnstructuredGrid") == 1 );
    REQUIRE( grid->GetNumberOfCells() == 1 );
    REQUIRE( grid->GetCellType(0) == VTK_HEXAHEDRON );
    REQUIRE( grid->HasAnyGhostCells() == false );

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

  }

  SECTION( "Check ghost cells" ) {

    short use_ghost_mask = 1;

    // Call CreateGrid function
    adaptor_creategrid(griddef.coords, griddef.npoints, griddef.cell_points, griddef.ncells,
                       griddef.ghost_mask, use_ghost_mask, mirror_periodic);

    vtkCPInputDataDescription * InputDescription = dataDescription->GetInputDescription(0);
    vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(InputDescription->GetGrid());

    REQUIRE( grid->HasAnyGhostCells() == true );

  }

  // Clean up test object and restore original one
  dataDescription->Delete();
  dataDescription = dataDescription_save;
  dataDescription_save = nullptr;

}

TEST_CASE( "Adaptor mirror_points with six cell periodic grid", "[coprocessor]" ) {

  TestGridDef griddef = NewSixCellPeriodicGrid();

  // New grid object
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Set up point coordinate array
  vtkSmartPointer<vtkDoubleArray> pointArray = vtkSmartPointer<vtkDoubleArray>::New();
  pointArray->SetNumberOfComponents(3);
  pointArray->SetNumberOfTuples(griddef.npoints);
  // Initialise with dummy coordinates, actual coords will be set later
  const double dummy_coords[] = {0.0, 0.0, 0.0};
  for (vtkIdType i = 0; i < griddef.npoints; i++) {
    pointArray->SetTuple(i, dummy_coords);
  }

  // Set up VTK points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(griddef.npoints);
  points->SetData(pointArray);

  // Hand over points to grid and define cells
  grid->SetPoints(points);
  grid->Allocate(static_cast<vtkIdType>(griddef.ncells));
  for(long cell = 0; cell < griddef.ncells; cell++) {
    vtkIdType thiscell[8] = {griddef.cell_points[cell*8  ], griddef.cell_points[cell*8+1], griddef.cell_points[cell*8+2], griddef.cell_points[cell*8+3],
                             griddef.cell_points[cell*8+4], griddef.cell_points[cell*8+5], griddef.cell_points[cell*8+6], griddef.cell_points[cell*8+7]};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, thiscell);
  }

  SECTION( "No mirror points for all xy > 0" ) {

    for (vtkIdType i = 0; i < griddef.npoints; i++) {
      pointArray->SetTuple(i, &griddef.coords[i*3]);
    }

    mirror_points(grid);

    REQUIRE( grid->GetNumberOfCells() == griddef.ncells );
    REQUIRE( grid->GetNumberOfPoints() == griddef.npoints );

  }

  SECTION( "Mirror points for xy shifted grid" ) {

    // Shift grid to set origin at grid center
    for (vtkIdType i = 0; i < griddef.npoints; i++) {
      const double newcoords[] = {griddef.coords[i*3]-1.5, griddef.coords[i*3+1]-1.5, griddef.coords[i*3+2]-0.5};
      pointArray->SetTuple(i, newcoords);
    }

    mirror_points(grid);

    // Same number of cells as before, but with 6 mirrored points
    REQUIRE( grid->GetNumberOfCells() == griddef.ncells );
    REQUIRE( grid->GetNumberOfPoints() == (griddef.npoints+6) );

    // Check that all cells have unit size
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); cellId++) {

      double cellBounds[6];
      grid->GetCellBounds(cellId, cellBounds);
      double cellDx = cellBounds[1] - cellBounds[0];
      double cellDy = cellBounds[3] - cellBounds[2];
      double cellDz = cellBounds[5] - cellBounds[4];

      REQUIRE ( cellDx == Approx(1.0) );
      REQUIRE ( cellDy == Approx(1.0) );
      REQUIRE ( cellDz == Approx(1.0) );

    }
  }
}

TEST_CASE( "Adaptor mirror_points with nine cell biperiodic grid", "[coprocessor]" ) {

  TestGridDef griddef = NewNineCellBiPeriodicGrid();

  // New grid object
  vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

  // Set up point coordinate array
  vtkSmartPointer<vtkDoubleArray> pointArray = vtkSmartPointer<vtkDoubleArray>::New();
  pointArray->SetNumberOfComponents(3);
  pointArray->SetNumberOfTuples(griddef.npoints);
  // Initialise with dummy coordinates, actual coords will be set later
  const double dummy_coords[] = {0.0, 0.0, 0.0};
  for (vtkIdType i = 0; i < griddef.npoints; i++) {
    pointArray->SetTuple(i, dummy_coords);
  }

  // Set up VTK points
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(griddef.npoints);
  points->SetData(pointArray);

  // Hand over points to grid and define cells
  grid->SetPoints(points);
  grid->Allocate(static_cast<vtkIdType>(griddef.ncells));
  for(long cell = 0; cell < griddef.ncells; cell++) {
    vtkIdType thiscell[8] = {griddef.cell_points[cell*8  ], griddef.cell_points[cell*8+1], griddef.cell_points[cell*8+2], griddef.cell_points[cell*8+3],
                             griddef.cell_points[cell*8+4], griddef.cell_points[cell*8+5], griddef.cell_points[cell*8+6], griddef.cell_points[cell*8+7]};
    grid->InsertNextCell(VTK_HEXAHEDRON, 8, thiscell);
  }

  SECTION( "No mirror points for all xy > 0" ) {

    for (vtkIdType i = 0; i < griddef.npoints; i++) {
      pointArray->SetTuple(i, &griddef.coords[i*3]);
    }

    mirror_points(grid);

    REQUIRE( grid->GetNumberOfCells() == griddef.ncells );
    REQUIRE( grid->GetNumberOfPoints() == griddef.npoints );

  }

  SECTION( "Mirror points for xy shifted grid" ) {

    // Shift grid to set origin at grid center
    for (vtkIdType i = 0; i < griddef.npoints; i++) {
      const double newcoords[] = {griddef.coords[i*3]-1.5, griddef.coords[i*3+1]-1.5, griddef.coords[i*3+2]-0.5};
      pointArray->SetTuple(i, newcoords);
    }

    mirror_points(grid);

    // Same number of cells as before, but with 14 mirrored points
    REQUIRE( grid->GetNumberOfCells() == griddef.ncells );
    REQUIRE( grid->GetNumberOfPoints() == (griddef.npoints+14) );

    // Check that all cells have unit size
    for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); cellId++) {

      double cellBounds[6];
      grid->GetCellBounds(cellId, cellBounds);
      double cellDx = cellBounds[1] - cellBounds[0];
      double cellDy = cellBounds[3] - cellBounds[2];
      double cellDz = cellBounds[5] - cellBounds[4];

      REQUIRE ( cellDx == Approx(1.0) );
      REQUIRE ( cellDy == Approx(1.0) );
      REQUIRE ( cellDz == Approx(1.0) );

    }
  }
}
