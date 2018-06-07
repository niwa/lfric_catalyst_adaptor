#include <unordered_map>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkCellType.h>
#include <vtkCellData.h>
#include <vtkUnsignedCharArray.h>

//
// Create Catalyst adaptor
//

// The purpose of the adaptor is to translate the simulation grid to
// a VTK grid, and to make the simulation data available that is stored
// in Fortran objects which are not directly accessible to C++ code.
// The easiest way to do this is to use simple buffers to move grid and
// data. The implementation uses VTK's smart pointers and stores the
// VTK grid and data to avoid scope issues and memory leaks.

// We need to access the dataDescription object defined in the coprocessor API
extern vtkCPDataDescription * dataDescription;

// This function mirrors points in periodic grids, where boundary cells reach
// to the opposite boundary. Works both for serial and parallel runs with
// grid partitioning. In the latter case, vertices (points) that are shared
// between partitions will be mirrored separately in each partition.
void mirror_points(vtkSmartPointer<vtkUnstructuredGrid> grid) {

  // Compute xy grid dimensions
  double gridBounds[6];
  grid->GetBounds(gridBounds);
  double gridDx = gridBounds[1] - gridBounds[0];
  double gridDy = gridBounds[3] - gridBounds[2];

  // Need to keep track of points that have already been duplicated,
  // to avoid degenerate points
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsX;
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsY;
  std::unordered_map<vtkIdType, vtkIdType> mirrorPointsXY;
  std::unordered_map<vtkIdType, vtkIdType>::const_iterator mirrorPointsIt;

  vtkPoints * gridPoints = grid->GetPoints();
  vtkSmartPointer<vtkIdList> oldCellPoints = vtkSmartPointer<vtkIdList>::New();

  // Search entire grid
  for (vtkIdType cellId = 0; cellId < grid->GetNumberOfCells(); cellId++) {

    // Compute xy cell dimensions
    double cellBounds[6];
    grid->GetCellBounds(cellId, cellBounds);
    double cellDx = cellBounds[1] - cellBounds[0];
    double cellDy = cellBounds[3] - cellBounds[2];

    // Find cells that span across the grid
    bool spanX = cellDx > 0.5*gridDx;
    bool spanY = cellDy > 0.5*gridDy;

    if (spanX or spanY) {

      grid->GetCellPoints(cellId, oldCellPoints);

      vtkIdType newCellPoints[8];

      // Check each cell vertex and mirror if needed
      for (vtkIdType pointIdIndex = 0; pointIdIndex < 8; pointIdIndex++) {

	vtkIdType thisPointId = oldCellPoints->GetId(pointIdIndex);
	double thisPointCoords[3];
	grid->GetPoint(thisPointId, thisPointCoords);

	// Mirror corner point
	if (spanX and spanY and thisPointCoords[0] < 0 and thisPointCoords[1] < 0) {
	  // Keep track of mirrored points to avoid degeneracy; insert a new point if
	  // no mirror point has been created yet
	  mirrorPointsIt = mirrorPointsXY.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsXY.end()) {
	    vtkIdType newPointId = gridPoints->InsertNextPoint(-thisPointCoords[0], -thisPointCoords[1], thisPointCoords[2]);
            mirrorPointsXY.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;   
	  }
          else {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// Mirror point on left domain boundary
	else if (spanX && thisPointCoords[0] < 0) {
          mirrorPointsIt = mirrorPointsX.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsX.end()) {
	    vtkIdType newPointId = gridPoints->InsertNextPoint(-thisPointCoords[0], thisPointCoords[1], thisPointCoords[2]);
            mirrorPointsX.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;   
	  }
          else {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// Mirror points on bottom domain boundary
	else if (spanY && thisPointCoords[1] < 0) {
          mirrorPointsIt = mirrorPointsY.find(thisPointId);
	  if (mirrorPointsIt == mirrorPointsY.end()) {
	    vtkIdType newPointId = gridPoints->InsertNextPoint(thisPointCoords[0], -thisPointCoords[1], thisPointCoords[2]);
            mirrorPointsY.insert({thisPointId, newPointId});
            newCellPoints[pointIdIndex] = newPointId;   
	  }
          else {
            newCellPoints[pointIdIndex] = mirrorPointsIt->second;
          }
	}
	// No mirror point needed
	else {
	  newCellPoints[pointIdIndex] = thisPointId;
	}

      }
      grid->ReplaceCell(cellId, 8, newCellPoints);
    }
  }
}

extern "C" {

  // Create a new VTK grid and register it with the coprocessor. Catalyst
  // supports several grids ("inputs"), but only one is used here.
  //
  // Funtion arguments:
  // point_coords: cartesian coordinates of each grid vertex
  // npoints: number of vertices in grid
  // cell_points: point indices for each cell
  // ncells: number of cells in grid
  // ghost_mask: mask array for ghost cells
  // use_ghost_mask: flag for including model ghost cells in VTK grid
  // mirror_periodic: flag for periodic grids which require replication of points
  void adaptor_creategrid(const double * point_coords, const long npoints,
                          const long * cell_points, const long ncells,
                          const short * ghost_mask, const short use_ghost_mask,
                          const short mirror_periodic) {

    if (!dataDescription) {
      vtkGenericWarningMacro("adaptor_creategrid: Unable to access dataDescription.");
      return;
    }
    // We only use one VTK grid and set of fields at this point
    if (dataDescription->GetNumberOfInputDescriptions() != 1) {
      vtkGenericWarningMacro("adaptor_creategrid: expected exactly 1 input description.");
      return;
    }
    if (npoints < 1) {
      vtkGenericWarningMacro("adaptor_creategrid: Invalid number provided for npoints:" << npoints);
      return;
    }
    if (ncells < 1) {
      vtkGenericWarningMacro("adaptor_creategrid: Invalid number provided for ncells:" << ncells);
      return;
    }

    // Create new grid
    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();

    //
    // Set grid point coordinates
    //

    // Copy coordinates, point_coords buffer will disappear
    vtkSmartPointer<vtkDoubleArray> pointArray = vtkSmartPointer<vtkDoubleArray>::New();
    pointArray->SetNumberOfComponents(3);
    pointArray->SetNumberOfTuples(npoints);
    for (vtkIdType i = 0; i < npoints; i++) {
      pointArray->SetTuple(i, &point_coords[i*3]);
    }

    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(npoints);
    points->SetData(pointArray);

    grid->SetPoints(points);

    //
    // Define grid cells
    //

    // We need memory for storing 8 points and cell type for each cell
    grid->Allocate(static_cast<vtkIdType>(ncells*9));

    for(size_t cell = 0; cell < ncells; cell++) {
      vtkIdType thiscell[8] = {cell_points[cell*8  ], cell_points[cell*8+1], cell_points[cell*8+2], cell_points[cell*8+3],
			       cell_points[cell*8+4], cell_points[cell*8+5], cell_points[cell*8+6], cell_points[cell*8+7]};
      grid->InsertNextCell(VTK_HEXAHEDRON, 8, thiscell);
    }

    //
    // Mark ghost cells if requested
    // It needs to be determined whether it is a good idea to transfer
    // existing ghost cells from the simulation to VTK
    //

    if (use_ghost_mask) {

      // Ghost cell array is not automatically allocated
      grid->AllocateCellGhostArray();

      // Ghost cells are treated as "duplicate cells" in VTK
      vtkUnsignedCharArray * ghosts = grid->GetCellGhostArray();
      for(vtkIdType cell = 0; cell < grid->GetNumberOfCells(); cell++) {
	if( ghost_mask[cell] == 1 ) {
	  ghosts->SetValue(cell, vtkDataSetAttributes::DUPLICATECELL);
	}
      }

    }

    // Mirror points in periodic grids, cells will otherwise span across
    // the entire extent of the grid
    if (mirror_periodic) {
      mirror_points(grid);
    }

    // Register grid with the coprocessor
    dataDescription->GetInputDescription(0)->SetGrid(grid);

  }

  // Copy field data from the simulation to VTK data structure
  // Funtion arguments:
  // fieldname: string that holds the name of the simulation field
  // fieldtype: point data (1) or cell data (2)
  // ncomponents: number of components in each data tuple
  // ntuples: the number of tuples in the field (number of points or cells)
  // fieldvalues: the field data itself
  void adaptor_copyfield(const char * fieldname, const int fieldtype, const int ncomponents,
			 const long ntuples, const double * fieldvalues) {

    if (!dataDescription) {
      vtkGenericWarningMacro("adaptor_copyfield: Unable to access dataDescription.");
      return;
    }
    if (dataDescription->GetNumberOfInputDescriptions() != 1) {
      vtkGenericWarningMacro("adaptor_copyfield: expected exactly 1 input description.");
      return;
    }

    // Try to get grid from coprocessor
    vtkCPInputDataDescription * InputDescription = dataDescription->GetInputDescription(0);
    vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(InputDescription->GetGrid());
    if (!grid) {
      vtkGenericWarningMacro("adaptor_copyfield: Unable to access VTK grid.");
      return;
    }

    // Sanity checking
    if (!strcmp(fieldname,"")) {
      vtkGenericWarningMacro("adaptor_copyfield: No name provided for fieldname.");
      return;
    }
    if (fieldtype != 1) {
      vtkGenericWarningMacro("adaptor_copyfield: Only cell data is supported at this point.");
      return;
    }
    if (grid->GetNumberOfCells() != ntuples) {
      vtkGenericWarningMacro("adaptor_copyfield: Number of tuples does not match VTK grid.");
      return;
    }

    // Check if field is needed and copy data
    if (InputDescription->IsFieldNeeded(fieldname)) {
      vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
      // Need to set number of components *before* number of tuples to ensure
      // that a sufficient amount of memory is allocated
      field->SetNumberOfComponents(ncomponents);
      field->SetNumberOfTuples(ntuples);
      field->SetName(fieldname);
      for (int j = 0; j < ncomponents; j++) {
	for (vtkIdType i = 0; i < ntuples; i++) {
	  field->SetComponent(i, j, fieldvalues[j*ntuples + i]);
	}
      }
      // Adding the array should automatically remove and deallocate arrays
      // that were previously added with the same name. So we shouldn't have
      // to worry about creating a memory leak here.
      grid->GetCellData()->AddArray(field);
    }

  }

}
