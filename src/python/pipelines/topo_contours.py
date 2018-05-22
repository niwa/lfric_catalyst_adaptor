import paraview.simple as pvs
from paraview import coprocessing as cp

#
# ParaView Catalyst pipeline for extracting a 2D spherical slice
# of the 3D mesh and overlaying it with coastlines and isocontours.
# The scene is rendered into an image and stored in png format.
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

# Values at which contours shall be produced
contour_values = [1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 15]

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:

      # Read topographic data
      topo = pvs.XMLPolyDataReader(FileName="ETOPO_10min_Ice.vtp")

      # Scale data to just under radius of the sphere
      toposcaled = pvs.Calculator(Input=topo)
      toposcaled.CoordinateResults = 1
      toposcaled.Function = '%i*coords' % (sphere_radius*0.99)

      # Define source of pipeline
      grid = coprocessor.CreateProducer( datadescription, "input" )

      ghosts = pvs.GhostCellsGenerator(Input=grid)
      ghosts.BuildIfRequired = 0
      ghosts.MinimumNumberOfGhostLevels = 1

      # Create a spherical slice
      slice = pvs.Slice(Input=ghosts)
      slice.SliceType = 'Sphere'
      slice.Triangulatetheslice = 0
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Radius = sphere_radius

      # Convert cell data to point data, which is required for good contour results
      #
      # CAUTION: THIS FILTER AVERAGES DATA FROM ALL CELLS SURROUNDING A POINT,
      #          WHICH REDUCES ACCURACY
      cell2point = pvs.CellDatatoPointData(Input=slice)
      cell2point.PassCellData = 0
      cell2point.PieceInvariant = 0

      # Create contours
      # Note that the "tube" filter can be used to highlight contours if needed.
      contours = pvs.Contour(Input=cell2point)
      contours.Isosurfaces=contour_values
      contours.ContourBy = ['POINTS', 'rho']
      contours.PointMergeMethod = 'Uniform Binning'

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
      coprocessor.RegisterView(renderView, filename='topo_contours_%t.png', freq=1,
                               fittoscreen=1, magnification=1, width=800, height=800,
                               cinema={})

      # Create colour transfer function for field
      LUT = pvs.GetColorTransferFunction('altitude')
      # Use Wikipedia LUT as provided on http://www.earthmodels.org/data-and-tools/color-tables
      LUT.RGBPoints = [	
        -11000, 0.141176470588235, 0.149019607843137,0.686274509803922,
        -5499.999, 0.219607843137255, 0.227450980392157,0.764705882352941,
        -5500, 0.219607843137255, 0.227450980392157,0.764705882352941,
        -2999.999, 0.274509803921569, 0.282352941176471,0.83921568627451,
        -3000, 0.274509803921569, 0.282352941176471,0.83921568627451,
        -1999.999, 0.317647058823529, 0.4,0.850980392156863,
        -2000, 0.317647058823529, 0.4,0.850980392156863,
        -749.999, 0.392156862745098, 0.505882352941176,0.874509803921569,
        -750, 0.392156862745098, 0.505882352941176,0.874509803921569,
        -69.999, 0.513725490196078, 0.631372549019608,0.901960784313726,
        -70, 0.513725490196078, 0.631372549019608,0.901960784313726,
        -19.999, 0.643137254901961, 0.752941176470588,0.941176470588235,
        -20, 0.643137254901961, 0.752941176470588,0.941176470588235,
        0.001, 0.666666666666667, 0.784313725490196,1,
        0, 0, 0.380392156862745,0.27843137254902,
        50.001, 0.0627450980392157, 0.47843137254902,0.184313725490196,
        50, 0.0627450980392157, 0.47843137254902,0.184313725490196,
        500.001, 0.909803921568627, 0.843137254901961,0.490196078431373,
        500, 0.909803921568627, 0.843137254901961,0.490196078431373,
        1200.001, 0.631372549019608, 0.262745098039216,0,
        1200, 0.631372549019608, 0.262745098039216,0,
        1700.001, 0.509803921568627, 0.117647058823529,0.117647058823529,
        1700, 0.509803921568627, 0.117647058823529,0.117647058823529,
        2800.001, 0.431372549019608, 0.431372549019608,0.431372549019608,
        2800, 0.431372549019608, 0.431372549019608,0.431372549019608,
        4000.001, 1, 1, 1,
        4000, 1, 1, 1,
        6000.001, 1, 1, 1
        ]

      LUT.ScalarRangeInitialized = 1.0

      # Show topo data
      sphereDisplay = pvs.Show(toposcaled, renderView)
      sphereDisplay.Representation = 'Surface'
      sphereDisplay.ColorArrayName = ['POINTS', 'altitude']
      sphereDisplay.LookupTable = LUT
      sphereDisplay.Opacity = 1.0

      # Show surface and colour by field value (which is cell data) using lookup table
      sphereDisplay = pvs.Show(contours, renderView)
      sphereDisplay.Representation = 'Surface'
      sphereDisplay.Opacity = 1.0

    return Pipeline()

  class CoProcessor(cp.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {'input': [1,1,1]}
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
