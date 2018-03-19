#!/usr/bin/python

import argparse
import vtk
from vtk.util import numpy_support as vtknps
import numpy as np

def xyz2llr(xyz):
    """
       Transforms an array of cartesian xyz coordinates
       into an array of longitude, latitude, and radius
       coordinates.
       Longitude range: -pi...pi
       Latitude range:  -pi/2...pi/2
    """
    xysq = xyz[:,0]**2 + xyz[:,1]**2
    llr = np.empty_like(xyz)
    
    # Longitude
    llr[:,0] = np.where(xyz[:,0] == 0,
                        0.5*np.pi*np.sign(xyz[:,1]),
                        np.arctan2(xyz[:,1], xyz[:,0]))
    # Latitude
    llr[:,1] = np.where(xysq == 0,
                        0.5*np.pi*np.sign(xyz[:,2]),
                        np.arctan2(xyz[:,2], np.sqrt(xysq)))
    # Radius - set to zero
    llr[:,2] = 0
    return llr

def llrprojection(map_projection_name):
    """
        Return a vtkGeoTransform object to transform a
        vtkPoints object with longitude-latitude (+radius)
        coordinates using a map projection selected by
        string "map_projection_name".
    """
    
    # Source points
    latlonproj = vtk.vtkGeoProjection()
    latlonproj.SetName('latlon')
    latlonproj.SetCentralMeridian(0)

    # Projected target points
    targetproj = vtk.vtkGeoProjection()
    targetproj.SetName(map_projection_name)
    targetproj.SetCentralMeridian(0)

    # Set up vtkGeoTransform object
    transform = vtk.vtkGeoTransform()
    transform.SetSourceProjection(latlonproj)
    transform.SetDestinationProjection(targetproj)
    return transform

def list_projections():
    """
        Print a formatted list of available map projections
        to stdout.
    """
    proj = vtk.vtkGeoProjection()
    print("Index   Short Name   Description")
    print("=================================")
    for idx in range(0, proj.GetNumberOfProjections()):
        print("%3i     %10s   %s" % (idx, proj.GetProjectionName(idx),
                                     proj.GetProjectionDescription(idx).split("\n")[0]))

def fix_crossing_cells(grid, verbose):
    """
        Find cells that span across the grid due to the
        "date line" where the sphere was cut, and fix them
        by replicating vertices ("points") on the other side.

        Note: this needs to be done in lat-lon coordinates,
        to make sure that cells are identified reliably.
    """
    if (verbose==True):
        print("Fixing cells that span across grid")

    # Set up cell locator
    cellloc = vtk.vtkCellLocator()
    cellloc.SetDataSet(grid)
    cellloc.BuildLocator()

    # Longitudes range between -pi...pi
    # Find cells along a meridional line on the last cube face,
    # just beyond longitude 3/2*pi. This will include both local cells
    # as well as cells than span over to the first cube face.
    grid.ComputeBounds()
    gridBounds = grid.GetBounds()
    # The last cube face is 3/2pi away; add 0.1 to make sure we are
    # within the cube face rather than at the edge
    searchLon = gridBounds[0] +1.6*np.pi
    lineStart = [searchLon, gridBounds[2], gridBounds[4]]
    lineEnd = [searchLon, gridBounds[3], gridBounds[5]]
    findTolerance = 0
    crossingCells = vtk.vtkIdList()
    cellloc.FindCellsAlongLine(lineStart, lineEnd, findTolerance, crossingCells)

    if (verbose == True):
        print("Grid boundaries (lon, lat, rad) = (%f:%f,%f:%f,%f:%f)" % gridBounds)
        print("Found", crossingCells.GetNumberOfIds(), "cells along line at lon=%f" % searchLon)

    gridPoints = grid.GetPoints()

    # Replicate points from the first cube face and replace cells
    # that span from the last to the first face
    replaceCellList = []
    cellCounter = 0
    for cellIdIndex in range(0, crossingCells.GetNumberOfIds()):

        CrossingCellId = crossingCells.GetId(cellIdIndex)
        pointIdList = grid.GetCell(CrossingCellId).GetPointIds()

        cellListed = False

        for pointIdindex in range(0, pointIdList.GetNumberOfIds()):

            pointId = pointIdList.GetId(pointIdindex)
            coords = grid.GetPoint(pointId)

            # Check if point lies on the first cube face, mirror at the meridian
            if ( coords[0] < -0.5*np.pi):
                newPointId = gridPoints.InsertNextPoint(-coords[0], coords[1], coords[2])
                replaceCellList.append([CrossingCellId, pointId, newPointId])
                cellListed = True

        # Check if this cell needed fixing, increment counter
        if (cellListed == True):
            cellCounter += 1

    # Redefine cells that need fixing using the additional points
    for cellId, oldPointId, newPointId in replaceCellList:
        grid.ReplaceCellPoint(cellId, oldPointId, newPointId)

    if (verbose == True):
        print("Replicated %i points" % len(replaceCellList))
        print("Replaced %i cells" % cellCounter)

def compute_projection(input_filename, output_filename, proj_name, verbose):

    #
    # Read poly data
    #

    # Check suffix - ".pvtp" is a parallel summary file, ".vtp" is a single polydata file
    if (input_filename.split(".")[-1] == "pvtp"):
        if (verbose == True):
            print("Reading input XML parallel poly data VTK file: %s" % input_filename)
        reader = vtk.vtkXMLPPolyDataReader()
    elif (input_filename.split(".")[-1] == "vtp"):
        if (verbose == True):
            print("Reading input XML poly data VTK file: %s" % input_filename)
        reader = vtk.vtkXMLPolyDataReader()
    else:
        print("ERROR - unknown file format suffix: %s" %
              input_filename.split(".")[-1])
        return

    reader.SetFileName(input_filename)
    reader.Update()
    polysphere = reader.GetOutput()

    # Remove ghost cells (if any)
    polysphere.RemoveGhostCells()

    if (verbose == True):
        print("Read %i vertices" % polysphere.GetNumberOfPoints())

    #
    # Transform point (vertex) coordinates to lon-lat-radius coordinates
    #

    points = polysphere.GetPoints()
    llrcoords = xyz2llr(vtknps.vtk_to_numpy(points.GetData()))
    points.SetData(vtknps.numpy_to_vtk(llrcoords))

    # Fix cells that span across the meridian
    # This is most easily done in lat-lon coordinates
    fix_crossing_cells(polysphere, verbose)

    #
    # Apply map projection
    #

    if (verbose == True):
        print("Applying %s projection to vertices" % proj_name)

    projpoints = vtk.vtkPoints()
    llrprojection(proj_name).TransformPoints(points, projpoints)
    polysphere.SetPoints(projpoints)

    if (output_filename.split(".")[-1] != "vtp"):
        print("WARNING - it is advisable to use suffix vtp for output data files!")

    if (verbose == True):
        print("Writing output XML poly data VTK file: %s" % output_filename)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetInputData(polysphere)
    writer.SetFileName(output_filename)
    writer.Update()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
                                     Applies a map projection to a VTK XML polydata file containing a sphere
                                     and writes the output into a VTK XML polydata file for further processing.
                                     Note that only vertex coordinates are projected, no transformation is
                                     applied to field data.
                                     
                                     Projections are computed using the PROJ library via VTK's vtkGeoProjection
                                     class. Available projections can be queried using the list-projections flag.
                                     """)
    parser.add_argument("--verbose", help="display additional information", action="store_true")
    parser.add_argument("--projname", help="short name of projection (default: moll)", default="moll", type=str)
    parser.add_argument("--list-projections", help="list all available projections and exit", action="store_true")
    parser.add_argument("input_file", help="VTK XML polydata file with a 2D sphere", type=str)
    parser.add_argument("output_file", help="VTK XML polydata file with map projection", type=str)
    args = parser.parse_args()

    if (args.list_projections == True):
        list_projections()
    else:
        compute_projection(args.input_file, args.output_file, args.projname, args.verbose)
