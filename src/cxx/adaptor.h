#ifndef ADAPTOR_H
#define ADAPTOR_H

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>


//
// Create Catalyst adaptor
//

// The purpose of the adaptor is to translate the simulation grid to
// a VTK grid, and to make the simulation data available that is stored
// in Fortran objects which are not directly accessible to C++ code.
// The easiest way to do this is to use simple buffers to move grid and
// data. The implementation uses VTK's smart pointers and stores the
// VTK grid and data to avoid scope issues and memory leaks.

// This function mirrors points in periodic grids, where boundary cells reach
// to the opposite boundary. Works both for serial and parallel runs with
// grid partitioning. In the latter case, vertices (points) that are shared
// between partitions will be mirrored separately in each partition.
void mirror_points(vtkSmartPointer<vtkUnstructuredGrid> grid);

extern "C" {

  // Create a new VTK grid and register it with the coprocessor. Catalyst
  // supports several grids ("inputs"), but only one is used here.
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
                          const short mirror_periodic);

  // Copy field data from the simulation to VTK data structure
  // Funtion arguments:
  // fieldname: string that holds the name of the simulation field
  // fieldtype: point data (1) or cell data (2)
  // ncomponents: number of components in each data tuple
  // ntuples: the number of tuples in the field (number of points or cells)
  // fieldvalues: the field data itself
  void adaptor_copyfield(const char * fieldname, const int fieldtype, const int ncomponents,
			 const long ntuples, const double * fieldvalues);

}

#endif
