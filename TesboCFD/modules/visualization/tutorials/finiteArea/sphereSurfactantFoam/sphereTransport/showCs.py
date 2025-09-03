
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.5.2

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=6
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.5.2

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.5.2

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1310, 756]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.StereoType = 0
      renderView1.CameraPosition = [-2.2544226375406193, 5.894774264129292, 2.2504803852287485]
      renderView1.CameraFocalPoint = [-0.17575377892067334, 0.05729450607871205, -0.2768849507248347]
      renderView1.CameraViewUp = [0.24075355565216316, -0.312162886467996, 0.9190168974250745]
      renderView1.CameraParallelScale = 0.9776975264012226
      renderView1.CameraParallelProjection = 1
      renderView1.Background = [0.32, 0.34, 0.43]

      # init the 'GridAxes3DActor' selected for 'AxesGrid'
      renderView1.AxesGrid.XTitleFontFile = ''
      renderView1.AxesGrid.YTitleFontFile = ''
      renderView1.AxesGrid.ZTitleFontFile = ''
      renderView1.AxesGrid.XLabelFontFile = ''
      renderView1.AxesGrid.YLabelFontFile = ''
      renderView1.AxesGrid.ZLabelFontFile = ''

      # register the view with coprocessor
      # and provide it with information such as the filename to use,
      # how frequently to write the images, etc.
      coprocessor.RegisterView(renderView1,
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1310, height=756, cinema={})
      renderView1.ViewTime = datadescription.GetTime()

      # ----------------------------------------------------------------
      # restore active view
      SetActiveView(renderView1)
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # setup the data processing pipelines
      # ----------------------------------------------------------------

      # create a new 'XML MultiBlock Data Reader'
      # create a producer from a simulation input
      area = coprocessor.CreateProducer(datadescription, 'area')

      # create a new 'Clip'
      clip1 = Clip(Input=area)
      clip1.ClipType = 'Plane'
      clip1.Scalars = ['CELLS', 'Cs']
      clip1.Value = 1.020149976015091

      # init the 'Plane' selected for 'ClipType'
      clip1.ClipType.Normal = [0.0, 0.0, 1.0]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from clip1
      clip1Display = Show(clip1, renderView1)

      # get color transfer function/color map for 'Cs'
      csLUT = GetColorTransferFunction('Cs')
      csLUT.RGBPoints = [0.7730422616004944, 0.231373, 0.298039, 0.752941, 1.020149976015091, 0.865003, 0.865003, 0.865003, 1.2672576904296875, 0.705882, 0.0156863, 0.14902]
      csLUT.ScalarRangeInitialized = 1.0

      # get opacity transfer function/opacity map for 'Cs'
      csPWF = GetOpacityTransferFunction('Cs')
      csPWF.Points = [0.7730422616004944, 0.0, 0.5, 0.0, 1.2672576904296875, 1.0, 0.5, 0.0]
      csPWF.ScalarRangeInitialized = 1

      # trace defaults for the display properties.
      clip1Display.Representation = 'Surface'
      clip1Display.ColorArrayName = ['CELLS', 'Cs']
      clip1Display.LookupTable = csLUT
      clip1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      clip1Display.SelectOrientationVectors = 'None'
      clip1Display.ScaleFactor = 0.2
      clip1Display.SelectScaleArray = 'None'
      clip1Display.GlyphType = 'Arrow'
      clip1Display.GlyphTableIndexArray = 'None'
      clip1Display.GaussianRadius = 0.01
      clip1Display.SetScaleArray = [None, '']
      clip1Display.ScaleTransferFunction = 'PiecewiseFunction'
      clip1Display.OpacityArray = [None, '']
      clip1Display.OpacityTransferFunction = 'PiecewiseFunction'
      clip1Display.DataAxesGrid = 'GridAxesRepresentation'
      clip1Display.SelectionCellLabelFontFile = ''
      clip1Display.SelectionPointLabelFontFile = ''
      clip1Display.PolarAxes = 'PolarAxesRepresentation'
      clip1Display.ScalarOpacityFunction = csPWF
      clip1Display.ScalarOpacityUnitDistance = 0.44814047465571655

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      clip1Display.DataAxesGrid.XTitleFontFile = ''
      clip1Display.DataAxesGrid.YTitleFontFile = ''
      clip1Display.DataAxesGrid.ZTitleFontFile = ''
      clip1Display.DataAxesGrid.XLabelFontFile = ''
      clip1Display.DataAxesGrid.YLabelFontFile = ''
      clip1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      clip1Display.PolarAxes.PolarAxisTitleFontFile = ''
      clip1Display.PolarAxes.PolarAxisLabelFontFile = ''
      clip1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      clip1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for csLUT in view renderView1
      csLUTColorBar = GetScalarBar(csLUT, renderView1)
      csLUTColorBar.Title = 'Cs'
      csLUTColorBar.ComponentTitle = ''
      csLUTColorBar.TitleFontFile = ''
      csLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      csLUTColorBar.Visibility = 1

      # show color legend
      clip1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(clip1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'area': [1, 1]}
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
    coprocessor.WriteImages(datadescription, rescale_lookuptable=rescale_lookuptable,
        image_quality=0, padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
