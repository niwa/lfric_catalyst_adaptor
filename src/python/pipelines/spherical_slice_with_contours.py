import paraview.simple as pvs
from paraview import coprocessing as cp

#
# ParaView Catalyst pipeline for extracting a 2D spherical slice
# of the 3D mesh, computing contour lines, and and storing
# everything  as VTK polygonal data for further analysis with
# ParaView or Python, along with data fields.
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
#

# NOTE: This currently only runs correctly in serial mode, as
#       neighbour/ghost cells are needed for the
#       CellDataToPointData filter.

#
# Pipeline parameters
#

# Radius of the sphere in [m]
sphere_radius = 6371230

# Values at which contours shall be produced
contour_values = [1.15, 1.2, 1.25, 1.3, 1.35]

# Switch for writing full model output
write_full_output = True

# Use static data range here to avoid MPI reduction call for
# establishing global range
dataRange = [1.2, 1.5]

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

      # Create a spherical slice
      slice = pvs.Slice(Input=grid)
      slice.SliceType = 'Sphere'
      slice.Triangulatetheslice = 0
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Radius = sphere_radius

      # Generate ghost cells - needed by CellDatatoPointData filter
      ghosts = pvs.GhostCellsGenerator(Input=grid)
      ghosts.BuildIfRequired = 0
      ghosts.MinimumNumberOfGhostLevels = 1

      # Convert cell data to point data, which is required for good contour results
      # Request "piece invariance" to ensure consistent values at
      # partition boundaries.
      #
      # CAUTION: THIS FILTER AVERAGES DATA FROM ALL CELLS SURROUNDING A POINT,
      #          WHICH REDUCES ACCURACY
      cell2point = pvs.CellDatatoPointData(Input=ghosts)

      # Create contours
      # Note that the "tube" filter can be used to highlight contours if needed.
      contours = pvs.Contour(Input=cell2point)
      contours.Isosurfaces=contour_values
      contours.ContourBy = ['POINTS', 'rho']
      contours.PointMergeMethod = 'Uniform Binning'

      # Create writers for slice and contour data and register them with the pipeline
      # Note that slice and contours generate separate datasets, so they need to be
      # written to separate files.
      sliceWriter = pvs.XMLPPolyDataWriter(Input=slice, DataMode="Appended",
                                           CompressorType="ZLib")
      coprocessor.RegisterWriter(sliceWriter, filename='spherical_slice_%t.pvtp', freq=1)

      # Create a new render view
      renderView = pvs.CreateView('RenderView')
      renderView.ViewSize = [1500, 768]
      renderView.AxesGrid = 'GridAxes3DActor'
      renderView.StereoType = 0
      renderView.CameraPosition = [0.0, 1.0, 0.3]
      renderView.CameraViewUp = [0.0, 0.0, 1.0]
      renderView.CameraParallelScale = 1.0
      renderView.Background = [0.32, 0.34, 0.43]
      renderView.ViewTime = datadescription.GetTime()

      # Register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView, filename='contours_%t.png', freq=1,
                               fittoscreen=1, magnification=1, width=1500, height=768,
                               cinema={})

      # Create colour transfer function for field
      LUT = pvs.GetColorTransferFunction('rho')
      LUT.RGBPoints = [dataRange[0], 0.23, 0.30, 0.75, 0.5*sum(dataRange), 0.87, 0.87, 0.87, dataRange[1], 0.71, 0.016, 0.15]
      LUT.ScalarRangeInitialized = 1.0

      # Show surface and colour by field value (which is cell data) using lookup table
      sphereDisplay = pvs.Show(slice, renderView)
      sphereDisplay.Representation = 'Surface'
      sphereDisplay.ColorArrayName = ['CELLS', 'rho']
      sphereDisplay.LookupTable = LUT

      # Show coastlines
      contourDisplay = pvs.Show(contours, renderView)
      contourDisplay.Representation = 'Surface'
      contourDisplay.ColorArrayName = [None, '']
      contourDisplay.OSPRayScaleArray = 'theta '
      contourDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      contourDisplay.SelectOrientationVectors = 'None'
      contourDisplay.ScaleFactor = 1193042.2418936265
      contourDisplay.SelectScaleArray = 'None'
      contourDisplay.GlyphType = 'Arrow'
      contourDisplay.GlyphTableIndexArray = 'None'
      contourDisplay.DataAxesGrid = 'GridAxesRepresentation'
      contourDisplay.PolarAxes = 'PolarAxesRepresentation'
      contourDisplay.GaussianRadius = 596521.1209468133
      contourDisplay.SetScaleArray = ['POINTS', 'theta ']
      contourDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      contourDisplay.OpacityArray = ['POINTS', 'theta ']
      contourDisplay.OpacityTransferFunction = 'PiecewiseFunction'

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
