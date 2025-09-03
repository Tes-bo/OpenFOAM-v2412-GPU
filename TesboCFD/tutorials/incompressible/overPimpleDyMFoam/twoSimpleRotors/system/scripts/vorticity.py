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
      renderView1.CameraPosition = [0.009999999776482582, 0.004999999888241291, 0.05232050690623429]
      renderView1.CameraFocalPoint = [0.009999999776482582, 0.004999999888241291, 0.004999999888241291]
      renderView1.CameraParallelScale = 0.01224744844016408
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
          filename='vorticity_%t.png', freq=1, fittoscreen=0, magnification=1, width=1091, height=766, cinema={})
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

      # create a new 'Stream Tracer'
      streamTracer1 = StreamTracer(Input=slice1,
          SeedType='High Resolution Line Source')
      streamTracer1.Vectors = ['POINTS', 'U']
      streamTracer1.MaximumStreamlineLength = 0.019999999552965164

      # init the 'High Resolution Line Source' selected for 'SeedType'
      streamTracer1.SeedType.Point1 = [0.0, 0.0, 0.005]
      streamTracer1.SeedType.Point2 = [0.02, 0.01, 0.005]

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from streamTracer1
      streamTracer1Display = Show(streamTracer1, renderView1)

      # get color transfer function/color map for 'Vorticity'
      vorticityLUT = GetColorTransferFunction('Vorticity')
      vorticityLUT.RGBPoints = [0.0, 0.229806, 0.298718, 0.753683, 37.5, 0.303869, 0.406535, 0.844959, 75.0, 0.383013, 0.509419, 0.917388, 112.5, 0.466667, 0.604563, 0.968155, 150.0, 0.552953, 0.688929, 0.995376, 187.5, 0.639176, 0.7596, 0.998151, 225.0, 0.722193, 0.813953, 0.976575, 262.5, 0.798692, 0.849786, 0.931689, 300.0, 0.865395, 0.86541, 0.865396, 337.5, 0.924128, 0.827385, 0.774508, 375.0, 0.958853, 0.769768, 0.678008, 412.5, 0.969954, 0.694267, 0.579375, 450.0, 0.958003, 0.602842, 0.481776, 487.50000000000006, 0.923945, 0.497309, 0.38797, 525.0, 0.869187, 0.378313, 0.300267, 562.5, 0.795632, 0.241284, 0.220526, 600.0, 0.705673, 0.0155562, 0.150233]
      vorticityLUT.ColorSpace = 'Lab'
      vorticityLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      streamTracer1Display.Representation = 'Surface'
      streamTracer1Display.ColorArrayName = ['POINTS', 'Vorticity']
      streamTracer1Display.LookupTable = vorticityLUT
      streamTracer1Display.OSPRayScaleArray = 'AngularVelocity'
      streamTracer1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      streamTracer1Display.SelectOrientationVectors = 'Normals'
      streamTracer1Display.ScaleFactor = 0.001999993808567524
      streamTracer1Display.SelectScaleArray = 'AngularVelocity'
      streamTracer1Display.GlyphType = 'Arrow'
      streamTracer1Display.GlyphTableIndexArray = 'AngularVelocity'
      streamTracer1Display.GaussianRadius = 9.99996904283762e-05
      streamTracer1Display.SetScaleArray = ['POINTS', 'AngularVelocity']
      streamTracer1Display.ScaleTransferFunction = 'PiecewiseFunction'
      streamTracer1Display.OpacityArray = ['POINTS', 'AngularVelocity']
      streamTracer1Display.OpacityTransferFunction = 'PiecewiseFunction'
      streamTracer1Display.DataAxesGrid = 'GridAxesRepresentation'
      streamTracer1Display.SelectionCellLabelFontFile = ''
      streamTracer1Display.SelectionPointLabelFontFile = ''
      streamTracer1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      streamTracer1Display.ScaleTransferFunction.Points = [-1.1626180405813291e-11, 0.0, 0.5, 0.0, 1.7840937690112886e-11, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      streamTracer1Display.OpacityTransferFunction.Points = [-1.1626180405813291e-11, 0.0, 0.5, 0.0, 1.7840937690112886e-11, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      streamTracer1Display.DataAxesGrid.XTitleFontFile = ''
      streamTracer1Display.DataAxesGrid.YTitleFontFile = ''
      streamTracer1Display.DataAxesGrid.ZTitleFontFile = ''
      streamTracer1Display.DataAxesGrid.XLabelFontFile = ''
      streamTracer1Display.DataAxesGrid.YLabelFontFile = ''
      streamTracer1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      streamTracer1Display.PolarAxes.PolarAxisTitleFontFile = ''
      streamTracer1Display.PolarAxes.PolarAxisLabelFontFile = ''
      streamTracer1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      streamTracer1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for vorticityLUT in view renderView1
      vorticityLUTColorBar = GetScalarBar(vorticityLUT, renderView1)
      vorticityLUTColorBar.Title = 'Vorticity'
      vorticityLUTColorBar.ComponentTitle = 'Magnitude'
      vorticityLUTColorBar.TitleFontFile = ''
      vorticityLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      vorticityLUTColorBar.Visibility = 1

      # show color legend
      streamTracer1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'Vorticity'
      vorticityPWF = GetOpacityTransferFunction('Vorticity')
      vorticityPWF.Points = [0.0, 0.0, 0.5, 0.0, 600.0, 1.0, 0.5, 0.0]
      vorticityPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(streamTracer1)
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
coprocessor.EnableLiveVisualization(True, 1)

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
