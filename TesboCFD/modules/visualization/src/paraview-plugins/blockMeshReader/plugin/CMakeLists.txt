#-----------------------------------------------------------------------------
# OpenFOAM blockMesh reader plugin for ParaView GUI

qt5_wrap_cpp(MOC_SRCS pqFoamBlockMeshControls.h)

paraview_plugin_add_property_widget(
    KIND        GROUP_WIDGET
    TYPE        "openfoam_blockreader_general_controls"
    CLASS_NAME  pqFoamBlockMeshControls
    INTERFACES  IFACES0
    SOURCES     IFACES0_SRCS)

set(BUILD_SHARED_LIBS_BACKUP ${BUILD_SHARED_LIBS})
set(BUILD_SHARED_LIBS ON)

paraview_add_plugin(ParaFoamBlockReader
    VERSION "${foam_api}"
    UI_INTERFACES
        ${IFACES0}
    SOURCES
        ${IFACES0_SRCS}
        ${MOC_SRCS}
        pqFoamBlockMeshControls.cxx
    MODULES blockMeshReader
    MODULE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/reader/vtk.module"
    SERVER_MANAGER_XML ParaFoamBlockReader.xml
)
set(BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS_BACKUP})

target_link_libraries(ParaFoamBlockReader
  PUBLIC
  ParaView::pqApplicationComponents
)

#-----------------------------------------------------------------------------
