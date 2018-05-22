import paraview.simple as pvs
from paraview import coprocessing as cp

#
# ParaView Catalyst pipeline for extracting a 2D spherical slice
# of the 3D mesh and showing velocity field magnitude and vectors.
# Creates a rendered image of the slice, and stores it as an image
# file in png format.
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
#

#
# Pipeline parameters
#

# Radius of the sphere in [m]
sphere_radius = 6371230

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      grid = coprocessor.CreateProducer(datadescription, 'input')

      # Normalise radius - simplifies creating glyphs
      normalise = pvs.Calculator(Input=grid)
      normalise.CoordinateResults = 1
      normalise.Function = 'coords/%i' % sphere_radius

      # Visualise velocity field using arrows
      glyph = pvs.Glyph(Input=normalise, GlyphType='Arrow')
      glyph.Scalars = ['POINTS', 'None']
      glyph.Vectors = ['CELLS', 'u']
      glyph.ScaleFactor = 0.2
      glyph.GlyphTransform = 'Transform2'
      glyph.GlyphType.TipResolution = 12
      glyph.GlyphType.TipRadius = 0.05
      glyph.GlyphType.ShaftRadius = 0.015

      # Create a new 'Render View'
      renderView = pvs.CreateView('RenderView')
      renderView.ViewSize = [1500, 768]
      renderView.AxesGrid = 'GridAxes3DActor'
      renderView.StereoType = 0
      renderView.CameraPosition = [-5, -2, 4]
      renderView.CameraViewUp = [0.5, 0.3, 0.8]
      renderView.CameraParallelScale = 1.7
      renderView.Background = [0.32, 0.34, 0.43]

      # Register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView, filename='velocity_field_%t.png', freq=1,
                               fittoscreen=1, magnification=1, width=1500, height=768,
                               cinema={})
      renderView.ViewTime = datadescription.GetTime()

      # Create colour transfer function for velocity field
      uLUT = pvs.GetColorTransferFunction('u')
      uLUT.RGBPoints = [1.7, 0.23, 0.30, 0.75, 20.9, 0.87, 0.87, 0.87, 40.0, 0.71, 0.016, 0.15]
      uLUT.ScalarRangeInitialized = 1.0

      # Show velocity field magnitude
      velocitymagDisplay = pvs.Show(normalise, renderView)
      velocitymagDisplay.Representation = 'Surface'
      velocitymagDisplay.ColorArrayName = ['CELLS', 'u']
      velocitymagDisplay.LookupTable = uLUT

      # Show colour legend
      uLUTColorBar = pvs.GetScalarBar(uLUT, renderView)
      uLUTColorBar.Title = 'u'
      uLUTColorBar.ComponentTitle = 'Magnitude'
      velocitymagDisplay.SetScalarBarVisibility(renderView, True)

      # Show velocity field glyphs
      glyphDisplay = pvs.Show(glyph, renderView)
      glyphDisplay.Representation = 'Surface'
      glyphDisplay.ColorArrayName = [None, '']

    return Pipeline()

  class CoProcessor(cp.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # These are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False)

# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False,
                            image_quality=0, padding_amount=0)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
