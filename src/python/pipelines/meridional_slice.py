import paraview.simple as pvs
from paraview import coprocessing as cp
import vtk
from vtk.util import numpy_support as vtknps
import numpy as np

#
# ParaView Catalyst pipeline for extracting a 2D meridional slice
# of the 3D mesh and storing it as VTK unstructured grid data for
# further analysis with ParaView or Python, along with data fields.
# Coordinate transformation is computed using NumPy, which should
# be reasonably fast.
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
#

#
# IMPORTANT NOTE: for low resolution cubed-sphere grids, the output
#                 of this pipeline may look warped at certain longitudes
#

#
# Pipeline parameters
#

# Longitude of the slice in degrees with the following convention
# with respect to Cartesian coordinates:
# -180 <= longitude < 0: westward
# longitude = 0: prime longitude (x axis)
# 0 < longitude < 180: eastward
longitude = 5

# Switch for writing full model output
write_full_output = True

# Switch for converting Cartesian coordinate to
# latitudes and radii
convert_to_lat_rad = True

# ----------------------------------------------------------------------

# Convert longitude to radians
lon_rad = longitude/180.0*np.pi

def xyz2latrad(xyz):
  xysq = xyz[:,0]**2 + xyz[:,1]**2
  llr = np.empty_like(xyz)

  # Latitude in degrees
  llr[:,0] = np.where(xysq == 0,
                      0.5*np.pi*np.sign(xyz[:,2]),
                      np.arctan2(xyz[:,2], np.sqrt(xysq)))
  llr[:,0] *= 180.0/np.pi
  # Longitude - set to zero
  llr[:,1] = 0
  # Radius - normalise
  llr[:,2] = np.sqrt(xysq + xyz[:,2]**2)
  radius_range = [np.min(llr[:,2]), np.max(llr[:,2])]
  llr[:,2] -= radius_range[0]
  llr[:,2] *= 90.0/(radius_range[1] - radius_range[0])

  return llr

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      # Define source of pipeline
      grid = coprocessor.CreateProducer( datadescription, "input" )

      if (write_full_output == True):
        fullWriter = pvs.XMLPUnstructuredGridWriter(Input=grid, DataMode="Appended",
                                                    CompressorType="ZLib")
        coprocessor.RegisterWriter(fullWriter, filename='full_output_%t.pvtu', freq=1)

      # Create a plane slice with the normal perpendicular
      # to the desired longitude
      slice = pvs.Slice(Input=grid)
      slice.SliceType = 'Plane'
      slice.Triangulatetheslice = 0
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Origin = [0,0,0]
      slice.SliceType.Normal = [-np.sin(lon_rad), np.cos(lon_rad), 0]

      # Clip the slice so that the hemisphere remains that
      # contains the desired longitude
      hemisphere = pvs.Clip(Input=slice)
      hemisphere.ClipType = 'Plane'
      hemisphere.ClipType.Origin = [0,0,0]
      hemisphere.ClipType.Normal = [np.cos(lon_rad), np.sin(lon_rad), 0]
      hemisphere.Crinkleclip = 0

      # We now switch from ParaView to plain VTK to transform
      # xyz to lon-rad coordinates if required. Note that this
      # could also be done in ParaView itself using the
      # Calculator filter - this may be the method of choice
      # when running the pipeline on a remote ParaView Server.

      # Run the pipeline first to instantiate all objects and
      # retrieve its vtkUnstructuredGrid output.
      # Note that this may not work when using a remote server.
      hemisphere.UpdatePipeline()
      data = hemisphere.GetClientSideObject().GetOutputDataObject(0)

      # Coordinate conversion for all vertices
      if (convert_to_lat_rad == True):
        points = data.GetPoints()
        newPoints = vtk.vtkPoints()
        latradcoords = xyz2latrad(vtknps.vtk_to_numpy(points.GetData()))
        newPoints.SetData(vtknps.numpy_to_vtk(latradcoords))
        data.SetPoints(newPoints)

      # Back to ParaView:
      # Create new writer and set its input to the VTK data object, then
      # register it with the coprocessor
      sliceWriter = pvs.XMLPUnstructuredGridWriter(DataMode="Appended", CompressorType="ZLib")
      sliceWriter.GetClientSideObject().SetInputDataObject(data)
      coprocessor.RegisterWriter(sliceWriter, filename='meridional_slice_%t.pvtu', freq=1)

    return Pipeline()

  class CoProcessor(cp.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {'input': [1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(False)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # Output all fields and meshes if forced
    if datadescription.GetForceOutput() == True:
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # Default implementation, uses output frequencies set in pipeline
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
