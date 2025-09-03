from paraview.simple import *
from paraview import coprocessing

fileNamePadding=6

# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      region = coprocessor.CreateProducer(datadescription, 'region')
      walls = coprocessor.CreateProducer(datadescription, 'walls')

      slice1 = Slice(Input=region)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]
      slice1.SliceType.Normal = [0, 0, 1]

      sliceWriter = servermanager.writers.XMLMultiBlockDataWriter(Input=slice1)
      coprocessor.RegisterWriter(sliceWriter, filename='slice_%t.vtm', freq=10, paddingamount=fileNamePadding)

      regionWriter = servermanager.writers.XMLMultiBlockDataWriter(Input=region)
      coprocessor.RegisterWriter(regionWriter, filename='region_%t.vtm', freq=100, paddingamount=fileNamePadding)

      wallsWriter = servermanager.writers.XMLMultiBlockDataWriter(Input=walls)
      coprocessor.RegisterWriter(wallsWriter, filename='walls_%t.vtm', freq=25, paddingamount=fileNamePadding)

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  freqs = {'region': [10, 100], 'walls': [25]}
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
coprocessor.EnableLiveVisualization(True)


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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
