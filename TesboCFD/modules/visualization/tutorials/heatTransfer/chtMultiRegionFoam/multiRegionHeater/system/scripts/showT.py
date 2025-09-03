
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
      renderView1.ViewSize = [1227, 756]
      renderView1.AxesGrid = 'GridAxes3DActor'
      renderView1.OrientationAxesVisibility = 0
      renderView1.CenterOfRotation = [0.0, -0.0159999993629754, 0.0]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [-0.22987670461473275, 0.027908517670887454, 0.37475316483990634]
      renderView1.CameraFocalPoint = [0.00145289821931723, -0.013255103910183037, 0.0005962967658004288]
      renderView1.CameraViewUp = [0.20896067334705018, 0.9776847628875804, 0.0216319534902209]
      renderView1.CameraParallelScale = 0.0945044135920057
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
          filename='image_%t.png', freq=1, fittoscreen=0, magnification=1, width=1227, height=756, cinema={})
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

      # create a new 'XML MultiBlock Data Reader'
      # create a producer from a simulation input
      walls = coprocessor.CreateProducer(datadescription, 'walls')

      # create a new 'Slice'
      slice1 = Slice(Input=region)
      slice1.SliceType = 'Plane'
      slice1.SliceOffsetValues = [0.0]

      # init the 'Plane' selected for 'SliceType'
      slice1.SliceType.Normal = [0.0, 0.0, 1.0]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from walls
      wallsDisplay = Show(walls, renderView1)

      # get color transfer function/color map for 'T'
      tLUT = GetColorTransferFunction('T')
      tLUT.RGBPoints = [300.0, 0.231373, 0.298039, 0.752941, 400.0, 0.865003, 0.865003, 0.865003, 500.0, 0.705882, 0.0156863, 0.14902]
      tLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      wallsDisplay.Representation = 'Surface'
      wallsDisplay.ColorArrayName = ['POINTS', 'T']
      wallsDisplay.LookupTable = tLUT
      wallsDisplay.Opacity = 0.4
      wallsDisplay.OSPRayScaleArray = 'T'
      wallsDisplay.OSPRayScaleFunction = 'PiecewiseFunction'
      wallsDisplay.SelectOrientationVectors = 'None'
      wallsDisplay.ScaleFactor = 0.0200000002980232
      wallsDisplay.SelectScaleArray = 'None'
      wallsDisplay.GlyphType = 'Arrow'
      wallsDisplay.GlyphTableIndexArray = 'None'
      wallsDisplay.GaussianRadius = 0.00100000001490116
      wallsDisplay.SetScaleArray = ['POINTS', 'T']
      wallsDisplay.ScaleTransferFunction = 'PiecewiseFunction'
      wallsDisplay.OpacityArray = ['POINTS', 'T']
      wallsDisplay.OpacityTransferFunction = 'PiecewiseFunction'
      wallsDisplay.DataAxesGrid = 'GridAxesRepresentation'
      wallsDisplay.SelectionCellLabelFontFile = ''
      wallsDisplay.SelectionPointLabelFontFile = ''
      wallsDisplay.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      wallsDisplay.ScaleTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 451.4297180175781, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      wallsDisplay.OpacityTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 451.4297180175781, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      wallsDisplay.DataAxesGrid.XTitleFontFile = ''
      wallsDisplay.DataAxesGrid.YTitleFontFile = ''
      wallsDisplay.DataAxesGrid.ZTitleFontFile = ''
      wallsDisplay.DataAxesGrid.XLabelFontFile = ''
      wallsDisplay.DataAxesGrid.YLabelFontFile = ''
      wallsDisplay.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      wallsDisplay.PolarAxes.PolarAxisTitleFontFile = ''
      wallsDisplay.PolarAxes.PolarAxisLabelFontFile = ''
      wallsDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
      wallsDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # show data from slice1
      slice1Display = Show(slice1, renderView1)

      # trace defaults for the display properties.
      slice1Display.Representation = 'Surface'
      slice1Display.ColorArrayName = ['POINTS', 'T']
      slice1Display.LookupTable = tLUT
      slice1Display.Opacity = 0.8
      slice1Display.OSPRayScaleArray = 'T'
      slice1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      slice1Display.SelectOrientationVectors = 'None'
      slice1Display.ScaleFactor = 0.0100000001490116
      slice1Display.SelectScaleArray = 'None'
      slice1Display.GlyphType = 'Arrow'
      slice1Display.GlyphTableIndexArray = 'None'
      slice1Display.GaussianRadius = 0.000500000007450581
      slice1Display.SetScaleArray = ['POINTS', 'T']
      slice1Display.ScaleTransferFunction = 'PiecewiseFunction'
      slice1Display.OpacityArray = ['POINTS', 'T']
      slice1Display.OpacityTransferFunction = 'PiecewiseFunction'
      slice1Display.DataAxesGrid = 'GridAxesRepresentation'
      slice1Display.SelectionCellLabelFontFile = ''
      slice1Display.SelectionPointLabelFontFile = ''
      slice1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      slice1Display.ScaleTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      slice1Display.OpacityTransferFunction.Points = [300.0, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]

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

      # get color legend/bar for tLUT in view renderView1
      tLUTColorBar = GetScalarBar(tLUT, renderView1)
      tLUTColorBar.Title = 'T'
      tLUTColorBar.ComponentTitle = ''
      tLUTColorBar.TitleFontFile = ''
      tLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      tLUTColorBar.Visibility = 1

      # show color legend
      wallsDisplay.SetScalarBarVisibility(renderView1, True)

      # show color legend
      slice1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'T'
      tPWF = GetOpacityTransferFunction('T')
      tPWF.Points = [300.0, 0.0, 0.5, 0.0, 500.0, 1.0, 0.5, 0.0]
      tPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(walls)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'walls': [1], 'region': [1, 1]}
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
