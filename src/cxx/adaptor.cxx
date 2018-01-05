#include <vtkCPAdaptorAPI.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkCellType.h>
#include <vtkCellData.h>

//
// Create Catalyst adaptor
//

// The purpose of the adaptor is to translate the simulation grid to
// a VTK grid, and to make the simulation data available to VTK.
// Fortran objects which are not directly accessible to C++ code.
// The easiest way to do this is to use simple buffers to
// move grid and data from to Catalyst. The implementation uses
// VTK's smart pointers and stores the VTK grid and data to avoid scope
// issues and memory leaks.

extern "C" {

  // Create a new VTK grid and register it with the coprocessor
  // Funtion arguments:
  // point_coords: cartesian coordinates of each grid vertex
  // npoints: number of vertices in grid
  // cell_points: point indices for each cell
  // ncells: number of cells in grid
  void adaptor_creategrid(const double * point_coords, const long npoints,
                          const long * cell_points, const long ncells) {

    if (!vtkCPAdaptorAPI::GetCoProcessorData()) {
      vtkGenericWarningMacro("adaptor_creategrid: Unable to access CoProcessorData.");
      return;
    }

    if(vtkCPAdaptorAPI::GetCoProcessorData()->GetNumberOfInputDescriptions() != 1) {
      vtkGenericWarningMacro("adaptor_creategrid: expected exactly 1 input description.");
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

    // Register grid with the coprocessor
    vtkCPAdaptorAPI::GetCoProcessorData()->GetInputDescription(0)->SetGrid(grid);

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

    if (!vtkCPAdaptorAPI::GetCoProcessorData()) {
      vtkGenericWarningMacro("adaptor_copyfield: Unable to access CoProcessorData.");
      return;
    }

    if(vtkCPAdaptorAPI::GetCoProcessorData()->GetNumberOfInputDescriptions() != 1) {
      vtkGenericWarningMacro("adaptor_copyfield: Expected exactly 1 input description.");
      return;
    }

    // Try to get grid from coprocessor
    vtkCPInputDataDescription * InputDescription = vtkCPAdaptorAPI::GetCoProcessorData()->GetInputDescription(0);
    vtkUnstructuredGrid * grid = vtkUnstructuredGrid::SafeDownCast(InputDescription->GetGrid());
    if (!grid) {
      vtkGenericWarningMacro("adaptor_copyfield: Unable to access VTK grid.");
      return;
    }

    // Sanity checking
    if (fieldtype != 1 | ncomponents != 1) {
      vtkGenericWarningMacro("adaptor_copyfield: Only scalar cell data is supported at this point.");
    }

    if (grid->GetNumberOfCells() != ntuples) {
      vtkGenericWarningMacro("adaptor_copyfield: Number of tuples does not match VTK grid.");
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
      grid->GetCellData()->AddArray(field);
    }
  }

}
