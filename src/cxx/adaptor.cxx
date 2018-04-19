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
  void adaptor_creategrid(const double * point_coords, const long npoints,
                          const long * cell_points, const long ncells,
                          const short * ghost_mask, const short use_ghost_mask) {

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
    if (fieldtype != 1 || ncomponents != 1) {
      vtkGenericWarningMacro("adaptor_copyfield: Only scalar cell data is supported at this point.");
      return;
    }
    if (grid->GetNumberOfCells() != ntuples) {
      vtkGenericWarningMacro("adaptor_copyfield: Number of tuples does not match VTK grid.");
      return;
    }

    // Check if field is needed and copy data
    if (InputDescription->IsFieldNeeded(fieldname)) {
      vtkSmartPointer<vtkDoubleArray> field = vtkSmartPointer<vtkDoubleArray>::New();
      field->SetNumberOfTuples(grid->GetNumberOfCells());
      field->SetNumberOfComponents(1);
      field->SetName(fieldname);
      for (vtkIdType i = 0; i < grid->GetNumberOfCells(); i++) {
	field->SetComponent(i, 0, fieldvalues[i]);
      }
      // Adding the array should automatically remove and deallocate arrays
      // that were previously added with the same name. So we shouldn't have
      // to worry about creating a memory leak here.
      grid->GetCellData()->AddArray(field);
    }

  }

}
