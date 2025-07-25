from paraview.simple import *
from paraview import coprocessing

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.5.0

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=4
rescale_lookuptable=False


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      # state file generated using paraview version 5.5.0

      # ----------------------------------------------------------------
      # setup views used in the visualization
      # ----------------------------------------------------------------

      # trace generated using paraview version 5.5.0

      #### disable automatic camera reset on 'Show'
      paraview.simple._DisableFirstRenderCameraReset()

      # Create a new 'Render View'
      renderView1 = CreateView('RenderView')
      renderView1.ViewSize = [1091, 766]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.CenterOfRotation = [0.009999999776482582, 0.004999999888241291, 0.004999999888241291]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [0.009999999776482582, 0.004999999888241291, 0.04819751509880177]
      renderView1.CameraFocalPoint = [0.009999999776482582, 0.004999999888241291, 0.004999999888241291]
      renderView1.CameraParallelScale = 0.011180339637598877
      renderView1.CameraParallelProjection = 1
      renderView1.Background = [0, 0, 0]

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
          filename='press_%t.png', freq=1, fittoscreen=0, magnification=1, width=1091, height=766, cinema={})
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
      region = coprocessor.CreateProducer(datadescription, 'region')

      # create a new 'Slice'
      slice1 = Slice(Input=region)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Origin = [0.01, 0.005, 0.005]
      slice1.SliceType.Normal = [0.0, 0.0, 1.0]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from slice1
      slice1Display = Show(slice1, renderView1)

      # get color transfer function/color map for 'p'
      pLUT = GetColorTransferFunction('p')
      pLUT.RGBPoints = [-0.2227432131767273, 0.231373, 0.298039, 0.752941, 0.0011433586478233337, 0.865003, 0.865003, 0.865003, 0.22502993047237396, 0.705882, 0.0156863, 0.14902]
      pLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      slice1Display.Representation = 'Surface'
      slice1Display.ColorArrayName = ['POINTS', 'p']
      slice1Display.LookupTable = pLUT
      slice1Display.OSPRayScaleArray = 'U'
      slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 0.0019999999552965165
      slice1Display.SelectScaleArray = 'None'
      slice1Display.GlyphType = 'Arrow'
      slice1Display.GlyphTableIndexArray = 'None'
      slice1Display.GaussianRadius = 9.999999776482583e-05
      slice1Display.SetScaleArray = ['POINTS', 'U']
      slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      slice1Display.OpacityArray = ['POINTS', 'U']
      slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
      slice1Display.DataAxesGrid = 'GridAxesRepresentation'
      slice1Display.SelectionCellLabelFontFile = ''
      slice1Display.SelectionPointLabelFontFile = ''
      slice1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      slice1Display.ScaleTransferFunction.Points = [-0.2436095029115677, 0.0, 0.5, 0.0, 0.2753259241580963, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      slice1Display.OpacityTransferFunction.Points = [-0.2436095029115677, 0.0, 0.5, 0.0, 0.2753259241580963, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      slice1Display.DataAxesGrid.XTitleFontFile = ''
      slice1Display.DataAxesGrid.YTitleFontFile = ''
      slice1Display.DataAxesGrid.ZTitleFontFile = ''
      slice1Display.DataAxesGrid.XLabelFontFile = ''
      slice1Display.DataAxesGrid.YLabelFontFile = ''
      slice1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      slice1Display.PolarAxes.PolarAxisTitleFontFile = ''
      slice1Display.PolarAxes.PolarAxisLabelFontFile = ''
      slice1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      slice1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for pLUT in view renderView1
      pLUTColorBar = GetScalarBar(pLUT, renderView1)
      pLUTColorBar.Title = 'p'
      pLUTColorBar.ComponentTitle = ''
      pLUTColorBar.TitleFontFile = ''
      pLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      pLUTColorBar.Visibility = 1

      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'p'
      pPWF = GetOpacityTransferFunction('p')
      pPWF.Points = [-0.2227432131767273, 0.0, 0.5, 0.0, 0.22502993047237396, 1.0, 0.5, 0.0]
      pPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(slice1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'region': [1, 1, 1]}
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
