
from paraview.simple import *
from paraview import coprocessing


#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.5.2

#--------------------------------------------------------------
# Global screenshot output options
imageFileNamePadding=0
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
      renderView1.CenterOfRotation = [0.5, 0.5, 0.5]
      renderView1.StereoType = 0
      renderView1.CameraPosition = [1.8716405831674128, 1.847832336175357, 3.241910834080083]
      renderView1.CameraFocalPoint = [0.6252880309830952, 0.5161911805422628, 0.4366497818793457]
      renderView1.CameraViewUp = [0.04375835431915484, 0.8948481435873961, -0.44422067302802376]
      renderView1.CameraParallelScale = 0.8660254037844386
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

      # create a new 'OpenFOAMReader'
      # create a producer from a simulation input
      cavityfoam = coprocessor.CreateProducer(datadescription, 'region')

      # create a new 'Contour'
      contour1 = Contour(Input=cavityfoam)
      contour1.ContourBy = ['POINTS', 'p_NOFIELD']
      contour1.Isosurfaces = [0.05]
      contour1.PointMergeMethod = 'Uniform Binning'

      # ----------------------------------------------------------------
      # setup the visualization in view 'renderView1'
      # ----------------------------------------------------------------

      # show data from contour1
      contour1Display = Show(contour1, renderView1)

      # get color transfer function/color map for 'U'
      uLUT = GetColorTransferFunction('U')
      uLUT.RGBPoints = [0.0012199414294188375, 0.231373, 0.298039, 0.752941, 0.5006099707147095, 0.865003, 0.865003, 0.865003, 1.0, 0.705882, 0.0156863, 0.14902]
      uLUT.ScalarRangeInitialized = 1.0

      # trace defaults for the display properties.
      contour1Display.Representation = 'Surface'
      contour1Display.ColorArrayName = ['POINTS', 'U']
      contour1Display.LookupTable = uLUT
      contour1Display.OSPRayScaleArray = 'Normals'
      contour1Display.OSPRayScaleFunction = 'PiecewiseFunction'
      contour1Display.SelectOrientationVectors = 'U'
      contour1Display.ScaleFactor = 0.1
      contour1Display.SelectScaleArray = 'p'
      contour1Display.GlyphType = 'Arrow'
      contour1Display.GlyphTableIndexArray = 'p'
      contour1Display.GaussianRadius = 0.005
      contour1Display.SetScaleArray = ['POINTS', 'Normals']
      contour1Display.ScaleTransferFunction = 'PiecewiseFunction'
      contour1Display.OpacityArray = ['POINTS', 'Normals']
      contour1Display.OpacityTransferFunction = 'PiecewiseFunction'
      contour1Display.DataAxesGrid = 'GridAxesRepresentation'
      contour1Display.SelectionCellLabelFontFile = ''
      contour1Display.SelectionPointLabelFontFile = ''
      contour1Display.PolarAxes = 'PolarAxesRepresentation'

      # init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
      contour1Display.ScaleTransferFunction.Points = [-0.9980009198188782, 0.0, 0.5, 0.0, -0.18620361387729645, 1.0, 0.5, 0.0]

      # init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
      contour1Display.OpacityTransferFunction.Points = [-0.9980009198188782, 0.0, 0.5, 0.0, -0.18620361387729645, 1.0, 0.5, 0.0]

      # init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
      contour1Display.DataAxesGrid.XTitleFontFile = ''
      contour1Display.DataAxesGrid.YTitleFontFile = ''
      contour1Display.DataAxesGrid.ZTitleFontFile = ''
      contour1Display.DataAxesGrid.XLabelFontFile = ''
      contour1Display.DataAxesGrid.YLabelFontFile = ''
      contour1Display.DataAxesGrid.ZLabelFontFile = ''

      # init the 'PolarAxesRepresentation' selected for 'PolarAxes'
      contour1Display.PolarAxes.PolarAxisTitleFontFile = ''
      contour1Display.PolarAxes.PolarAxisLabelFontFile = ''
      contour1Display.PolarAxes.LastRadialAxisTextFontFile = ''
      contour1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''

      # setup the color legend parameters for each legend in this view

      # get color legend/bar for uLUT in view renderView1
      uLUTColorBar = GetScalarBar(uLUT, renderView1)
      uLUTColorBar.Title = 'U'
      uLUTColorBar.ComponentTitle = 'Magnitude'
      uLUTColorBar.TitleFontFile = ''
      uLUTColorBar.LabelFontFile = ''

      # set color bar visibility
      uLUTColorBar.Visibility = 1

      # show color legend
      contour1Display.SetScalarBarVisibility(renderView1, True)

      # ----------------------------------------------------------------
      # setup color maps and opacity mapes used in the visualization
      # note: the Get..() functions create a new object, if needed
      # ----------------------------------------------------------------

      # get opacity transfer function/opacity map for 'U'
      uPWF = GetOpacityTransferFunction('U')
      uPWF.Points = [0.0012199414294188375, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
      uPWF.ScalarRangeInitialized = 1

      # ----------------------------------------------------------------
      # finally, restore active source
      SetActiveSource(contour1)
      # ----------------------------------------------------------------
    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  coprocessor = CoProcessor()
  # these are the frequencies at which the coprocessor updates.
  freqs = {'region': [1, 1]}
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
