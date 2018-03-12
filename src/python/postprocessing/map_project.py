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
       Longitude range: 0...2pi
       Latitude range:  -pi/2...pi/2
    """
    xysq = xyz[:,0]**2 + xyz[:,1]**2
    llr = np.empty_like(xyz)
    
    # Longitude
    llr[:,0] = np.where(xyz[:,0] == 0,
                        0.5*np.pi*np.sign(xyz[:,1]),
                        np.arctan2(xyz[:,1], xyz[:,0]))
    llr[:,0] += np.where(llr[:,0] < 0, 2*np.pi, 0)
    # Latitude
    llr[:,1] = np.where(xysq == 0,
                        0.5*np.pi*np.sign(xyz[:,2]),
                        np.arctan2(xyz[:,2], np.sqrt(xysq)))
    # Radius
    llr[:,2] = np.sqrt(xysq + xyz[:,2]**2)
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

    if (verbose == True):
        print("Read %i vertices" % polysphere.GetNumberOfPoints())

    #
    # Transform point (vertex) coordinates to lon-lat-radius coordinates
    #

    points = polysphere.GetPoints()
    llrcoords = xyz2llr(vtknps.vtk_to_numpy(points.GetData()))
    points.SetData(vtknps.numpy_to_vtk(llrcoords))

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
