import paraview.simple as pvs
from paraview import coprocessing as cp

#
# ParaView Catalyst pipeline for extracting a 2D spherical slice
# of the 3D mesh and storing it as VTK polygonal data for further
# analysis with ParaView or Python, along with data fields. This
# version overlays coastlines, creates a rendered image of the
# slice, and stores it as an image file in png format.
#
# Coastlines data set can be downloaded here:
# http://www.earthmodels.org/data-and-tools/coastlines/Coastlines_Los_Alamos.vtp
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
#

#
# IMPORTANT NOTE: this pipeline requires rendering components
#                 that may require a full build of ParaView
#

#
# Pipeline parameters
#

# Radius of the sphere in [m]
sphere_radius = 6371230

# Switch for writing full model output
write_full_output = True

# Name of field used for colouring the rendered image
fieldname = 'theta'

# Use static data range here to avoid MPI reduction call for
# establishing global range
dataRange = [2.7e2, 2.8e2]

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      # Read coastlines
      clines = pvs.XMLPolyDataReader(FileName="Coastlines_Los_Alamos.vtp")

      # Scale data to radius of the sphere
      clinesscaled = pvs.Calculator(Input=clines)
      clinesscaled.CoordinateResults = 1
      clinesscaled.Function = '%i*coords' % sphere_radius

      # Define source of pipeline
      grid = coprocessor.CreateProducer( datadescription, "input" )

      if (write_full_output == True):
        fullWriter = pvs.XMLPUnstructuredGridWriter(Input=grid, DataMode="Appended",
                                                    CompressorType="ZLib")
        coprocessor.RegisterWriter(fullWriter, filename='full_output_%t.pvtu', freq=1)

      # Create a spherical slice
      slice = pvs.Slice(Input=grid)
      slice.SliceType = 'Sphere'
      slice.Triangulatetheslice = 0
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Radius = sphere_radius

      # Create writer for this data and register it with the pipeline
      sliceWriter = pvs.XMLPPolyDataWriter(Input=slice, DataMode="Appended",
                                           CompressorType="ZLib")
      coprocessor.RegisterWriter(sliceWriter, filename='spherical_slice_%t.pvtp', freq=1)

      # Create a new render view
      renderView = pvs.CreateView('RenderView')
      renderView.ViewSize = [1500, 768]
      renderView.AxesGrid = 'GridAxes3DActor'
      renderView.StereoType = 0
      renderView.CameraPosition = [-5, -2, 4]
      renderView.CameraViewUp = [0.5, 0.3, 0.8]
      renderView.CameraParallelScale = 1.7
      renderView.Background = [0.32, 0.34, 0.43]
      renderView.ViewTime = datadescription.GetTime()

      # Register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView, filename='spherical_slice_%t.png', freq=1,
                               fittoscreen=1, magnification=1, width=800, height=800,
                               cinema={})

      # Create colour transfer function for field
      LUT = pvs.GetColorTransferFunction(fieldname)
      LUT.RGBPoints = [dataRange[0], 0.23, 0.30, 0.75, 0.5*sum(dataRange), 0.87, 0.87, 0.87, dataRange[1], 0.71, 0.016, 0.15]
      LUT.ScalarRangeInitialized = 1.0

      # Show surface and colour by field value (which is cell data) using lookup table
      sphereDisplay = pvs.Show(slice, renderView)
      sphereDisplay.Representation = 'Surface'
      sphereDisplay.ColorArrayName = ['CELLS', fieldname]
      sphereDisplay.LookupTable = LUT
      sphereDisplay.Opacity = 1.0

      # Show coastlines
      sphereDisplay = pvs.Show(clinesscaled, renderView)
      sphereDisplay.Representation = 'Surface'
      sphereDisplay.Opacity = 1.0

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
