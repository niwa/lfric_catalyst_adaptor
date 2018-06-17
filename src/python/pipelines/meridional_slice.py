import paraview.simple as pvs
from paraview import coprocessing as cp
import math

#
# ParaView Catalyst pipeline for extracting a 2D meridional slice
# of the 3D mesh and storing it as VTK unstructured grid data for
# further analysis with ParaView or Python, along with data fields.
#
# A list of available filters and writers can be found here:
# https://www.paraview.org/ParaView/Doc/Nightly/www/py-doc/\
# paraview.servermanager_proxies.html
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

      # Determine grid bounds
      gridBounds = grid.GetClientSideObject().GetOutputDataObject(0).GetBounds()

      # Plane slice at the desired longitude
      slice = pvs.Slice(Input=grid)
      slice.SliceType = 'Plane'
      slice.Triangulatetheslice = 0
      slice.SliceOffsetValues = [0.0]
      slice.SliceType.Origin = [longitude/180.0*math.pi, 0, gridBounds[4]]
      slice.SliceType.Normal = [1, 0, 0]

      # Create new writer
      sliceWriter = pvs.XMLPPolyDataWriter(Input=slice, DataMode="Appended", CompressorType="ZLib")
      coprocessor.RegisterWriter(sliceWriter, filename='meridional_slice_%t.pvtp', freq=1)

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
