import paraview.simple as pvs
from paraview import coprocessing as cp

#
# ParaView Catalyst pipeline for visualising global winds in 3D using
# lon lat rad coordinates. Creates a rendered image of the slice, and
# stores it as an image file in png format.
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
#

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      grid = coprocessor.CreateProducer(datadescription, 'input')

      # Simulation domain outline
      outline = pvs.Outline(Input=grid)

      # Horizontal slice at the bottom
      slice = pvs.Slice(Input=grid)
      slice.SliceType = 'Plane'
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Origin = [0.0, 0.0, 101.94]
      slice.SliceType.Normal = [0.0, 0.0, 1.0]
      slice.Triangulatetheslice = False

      # Glyphs for representing velocity field
      glyph = pvs.Glyph(Input=grid, GlyphType='Arrow')
      glyph.Vectors = ['CELLS', 'u']
      glyph.ScaleMode = 'vector'
      glyph.ScaleFactor = 0.01
      glyph.GlyphMode = 'Every Nth Point'
      glyph.Stride = 200

      # Create a new render view
      renderView = pvs.CreateView('RenderView')
      renderView.ViewSize = [800, 400]
      renderView.InteractionMode = '2D'
      renderView.AxesGrid = 'GridAxes3DActor'
      renderView.CenterOfRotation = [0.18, 0.0, 102]
      renderView.StereoType = 0
      renderView.CameraPosition = [-3.4, -6.8, 107]
      renderView.CameraFocalPoint = [-0.27, -0.41, 102]
      renderView.CameraViewUp = [0.057, 0.49, 0.87]
      renderView.CameraParallelScale = 1.0
      renderView.Background = [0.32, 0.34, 0.43]

      # Register the view with coprocessor
      coprocessor.RegisterView(renderView, filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=800, height=400, cinema={})
      renderView.ViewTime = datadescription.GetTime()

      # Get color transfer function/color map for field rho
      rhoLUT = pvs.GetColorTransferFunction('rho')
      rhoLUT.RGBPoints = [1.17, 0.231, 0.298, 0.752, 1.33, 0.865, 0.865, 0.865, 1.49, 0.706, 0.0157, 0.149]
      rhoLUT.ScalarRangeInitialized = 1.0

      # Show slice
      sliceDisplay = pvs.Show(slice, renderView)
      sliceDisplay.Representation = 'Surface With Edges'
      sliceDisplay.ColorArrayName = ['CELLS', 'rho']
      sliceDisplay.LookupTable = rhoLUT
      sliceDisplay.ScaleFactor = 0.628
      sliceDisplay.SelectScaleArray = 'None'
      sliceDisplay.GlyphType = 'Arrow'
      sliceDisplay.GlyphTableIndexArray = 'None'
      sliceDisplay.DataAxesGrid = 'GridAxesRepresentation'
      sliceDisplay.PolarAxes = 'PolarAxesRepresentation'
      sliceDisplay.GaussianRadius = 0.314
      sliceDisplay.SetScaleArray = [None, '']
      sliceDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      sliceDisplay.OpacityArray = [None, '']
      sliceDisplay.OpacityTransferFunction = 'PiecewiseFunction'

      # Show color legend
      sliceDisplay.SetScalarBarVisibility(renderView, True)

      # Show glyph
      glyphDisplay = pvs.Show(glyph, renderView)

      # Show outline
      outlineDisplay = pvs.Show(outline, renderView)

      # Get color legend/bar for rhoLUT in view renderView
      rhoLUTColorBar = pvs.GetScalarBar(rhoLUT, renderView)
      rhoLUTColorBar.WindowLocation = 'LowerRightCorner'
      rhoLUTColorBar.Title = 'rho'
      rhoLUTColorBar.ComponentTitle = ''

    return Pipeline()

  class CoProcessor(cp.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # These are the frequencies at which the coprocessor updates.
  freqs = {'input': [1, 1, 1, 1]}
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
coprocessor.EnableLiveVisualization(False, 1)

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
