#cmake_minimum_required(VERSION 2.6.3)

#FOREACH(policy CMP0012 CMP0013 CMP0014 CMP0017)
#  IF(POLICY ${policy})
#    CMAKE_POLICY(SET ${policy} NEW)
#  ENDIF()
#ENDFOREACH()

project(VisItBridgePlugin)

include(${CMAKE_CURRENT_SOURCE_DIR}/CMake/VisItCommon.cmake)

#-----------------------------------------------------------------------------
# Create visit-config.h
#-----------------------------------------------------------------------------

if(WIN32)
  add_definitions(-DWIN32)
endif()

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/include/visit-cmake.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/visit-config.h @ONLY IMMEDIATE)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/common/misc/ConfigureInfo.C.in
               ${CMAKE_CURRENT_BINARY_DIR}/common/misc/ConfigureInfo.C @ONLY IMMEDIATE)


add_subdirectory(common)
add_subdirectory(visit_vtk)
add_subdirectory(avt)

if (VISIT_PARALLEL)
  set(visit_par_libs
    avtdatabase_par
    avtpipeline_par)
else ()
  set(visit_par_libs
    avtdatabase_ser
    avtpipeline_ser)
endif ()

add_library(VisItLib INTERFACE)
target_link_libraries(VisItLib
  INTERFACE
    visitcommon
    ${visit_par_libs}
    visit_vtk)
target_include_directories(VisItLib
  INTERFACE
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
add_library(ParaView::VisItLib ALIAS VisItLib)
VISIT_INSTALL_TARGETS(VisItLib)
