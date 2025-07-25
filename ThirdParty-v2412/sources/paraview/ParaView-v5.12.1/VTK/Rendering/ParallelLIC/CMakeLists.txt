set(classes
  vtkMPIPixelTT
  vtkParallelTimer
  vtkPLineIntegralConvolution2D
  vtkPPainterCommunicator
  vtkPPixelTransfer
  vtkPSurfaceLICComposite
  vtkPSurfaceLICInterface)
set(headers
  vtkMPIPixelView.h)

vtk_object_factory_declare(
  BASE vtkSurfaceLICInterface
  OVERRIDE vtkPSurfaceLICInterface)
vtk_object_factory_declare(
  BASE vtkLineIntegralConvolution2D
  OVERRIDE vtkPLineIntegralConvolution2D)
vtk_object_factory_declare(
  BASE vtkSurfaceLICComposite
  OVERRIDE vtkPSurfaceLICComposite)

vtk_object_factory_configure(
  SOURCE_FILE vtk_object_factory_source
  HEADER_FILE vtk_object_factory_header
  EXPORT_MACRO "VTKRENDERINGPARALLELLIC_EXPORT")

set(shader_files
  vtkPSurfaceLICComposite_CompFS.glsl)

set(sources)
set(private_headers)
foreach(shader_file IN LISTS shader_files)
  vtk_encode_string(
    INPUT         "${shader_file}"
    HEADER_OUTPUT header
    SOURCE_OUTPUT source)
  list(APPEND sources ${source})
  list(APPEND private_headers ${header})
endforeach()

vtk_module_add_module(VTK::RenderingParallelLIC
  CLASSES ${classes}
  SOURCES ${sources} ${vtk_object_factory_source}
  HEADERS ${headers}
  PRIVATE_HEADERS ${private_headers} ${vtk_object_factory_header})

option(VTK_RENDERINGPARALLELLIC_SURFACELICPAINTER_TIMER "enable parallel timers for the surface lic painter" OFF)
mark_as_advanced(VTK_RENDERINGPARALLELLIC_SURFACELICPAINTER_TIMER)
vtk_deprecated_setting(convolution_default
  VTK_RENDERINGPARALLELLIC_LINEINTEGRALCONVOLUTION2D_TIMER
  VTK_RENDERINGPARALLELLIC_LINEINTEGRALCONVLOLUTION2D_TIMER
  OFF)
option(VTK_RENDERINGPARALLELLIC_LINEINTEGRALCONVOLUTION2D_TIMER "enable parallel timers for the 2d line integral convolution" "${convolution_default}")
mark_as_advanced(VTK_RENDERINGPARALLELLIC_LINEINTEGRALCONVOLUTION2D_TIMER)

if (VTK_RENDERINGPARALLELLIC_SURFACELICPAINTER_TIMER)
  vtk_module_definitions(VTK::RenderingParallelLIC
    PRIVATE
      vtkSurfaceLICPainterTIME)
endif()
if (VTK_RENDERINGPARALLELLIC_LINEINTEGRALCONVOLUTION2D_TIMER)
  vtk_module_definitions(VTK::RenderingParallelLIC
    PRIVATE
      vtkLineIntegralConvolution2DTIME)
endif()
vtk_add_test_mangling(VTK::RenderingParallelLIC)
