#get same cmake config as visitlib itself does
include(${CMAKE_CURRENT_SOURCE_DIR}/../Library/VisItLib/CMake/VisItCommon.cmake)

macro(VISIT_INCLUDE_READER_SOURCE name state)
option(VISIT_BUILD_READER_${name} "Build VisIt Reader ${name}" ${state})
if (VISIT_BUILD_READER_${name})
  add_subdirectory(readers/${name})
endif()
endmacro(VISIT_INCLUDE_READER_SOURCE name)

#these are used by visit_include_reader
set(VISIT_DATABASE_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})
set(VISIT_SERVER_SOURCES CACHE INTERNAL "vtk classes to wrap for client server" FORCE)
set(VISIT_DB_SOURCES "" CACHE INTERNAL "visit sources for readers" FORCE)
set(VISIT_DB_INC_DIRS CACHE INTERNAL "include directories" FORCE)

set(DEFAULT_BRIDGE_READERS
  ANSYS
  AUXFile
  BOV
  CEAucd
  Chombo
  CMAT
  CTRL
  Curve2D
  DDCMD
  Dyna3D
  ExtrudedVol
  FLASH
  Fluent
  Gadget
  GGCM
  GTC
  GULP
  H5Nimrod
  Image
  LAMMPS
  Lines
  M3D
  M3DC1
  MFIX
  MFIXCDF
  Miranda
  MM5
  NASTRAN
  NETCDF
  Nek5000
  OpenFOAM
  OVERFLOW
  paraDIS
  paraDIS_tecplot
  PATRAN
  PFLOTRAN
  Pixie
  PLOT3D
  Point3D
  ProteinDataBank
  RAW
  SAMRAI
  SAR
  SAS
  Spheral
  Tecplot
  Tetrad
  TFT
  TSurf
  UNIC
  VASP
  Velodyne
  Vs
  Xmdv
  XYZ
  )
if(NOT WIN32)
  list(APPEND DEFAULT_BRIDGE_READERS
    Claw
    )
endif(NOT WIN32)

#add generic vtkreader interface parent classes
add_subdirectory(AvtAlgorithms)

list(INSERT CMAKE_MODULE_PATH 0
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
set(_VisItReader_output_dir
  "${CMAKE_CURRENT_BINARY_DIR}")
include(VisItReader)

#add the default readers to the library
foreach(db ${DEFAULT_BRIDGE_READERS})
  add_subdirectory(readers/${db})
endforeach(db)

add_definitions("-DRC_CPP_VISIT_BUILD")

# set of libraries xmls for readers enabled.
set(visit_libraries)
set(visit_readers_xmls
  visit_readers.xml)

#-----------------------------------------------------------------------------
visit_include_reader_source(Silo OFF)
if(VISIT_BUILD_READER_Silo)
  vtk_module_find_package(PACKAGE SILO)
  list(APPEND visit_libraries
    SILO::SILO)
  list(APPEND visit_readers_xmls
    visit_silo_reader.xml)
endif()

#-----------------------------------------------------------------------------
if(VISIT_BUILD_READER_Boxlib3D)
  vtk_module_find_package(PACKAGE Boxlib)
  list(APPEND visit_libraries
    Boxlib::Boxlib)
  list(APPEND visit_readers_xmls
    visit_boxlib3d_reader.xml)
endif()
visit_include_reader_source(Boxlib3D OFF)

#-----------------------------------------------------------------------------
visit_include_reader_source(Mili OFF)
if(VISIT_BUILD_READER_Mili)
  vtk_module_find_package(PACKAGE Mili)
  vtk_module_find_package(PACKAGE RapidJSON)
  list(APPEND visit_libraries
    Mili::Mili
    rapidjson)
  list(APPEND visit_readers_xmls
    visit_mili_reader.xml)
endif()

set(sources
  AvtAlgorithms/vtkAvtFileFormatAlgorithm.cxx
  AvtAlgorithms/vtkAvtMTMDFileFormatAlgorithm.cxx
  AvtAlgorithms/vtkAvtMTSDFileFormatAlgorithm.cxx
  AvtAlgorithms/vtkAvtSTMDFileFormatAlgorithm.cxx
  AvtAlgorithms/vtkAvtSTSDFileFormatAlgorithm.cxx)

set(headers
  AvtAlgorithms/vtkAvtFileFormatAlgorithm.h
  AvtAlgorithms/vtkAvtMTMDFileFormatAlgorithm.h
  AvtAlgorithms/vtkAvtMTSDFileFormatAlgorithm.h
  AvtAlgorithms/vtkAvtSTMDFileFormatAlgorithm.h
  AvtAlgorithms/vtkAvtSTSDFileFormatAlgorithm.h)

set(private_headers
  AvtAlgorithms/vtkFileCloseSemiProperty.h)

get_property(visit_server_sources GLOBAL PROPERTY _visit_server_sources)
get_property(visit_db_sources GLOBAL PROPERTY _visit_db_sources)
get_property(visit_db_headers GLOBAL PROPERTY _visit_db_headers)
get_property(visit_include_dirs GLOBAL PROPERTY _visit_db_include_dirs)

vtk_module_add_module(ParaView::IOVisItBridge
  SPDX_SKIP_REGEX ".*"
  SOURCES ${visit_server_sources} ${visit_db_sources} ${sources}
  HEADERS ${headers} ${visit_db_headers}
  PRIVATE_HEADERS ${private_headers})
vtk_module_include(ParaView::IOVisItBridge
  PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/AvtAlgorithms>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
  PRIVATE
    ${visit_include_dirs})

paraview_add_server_manager_xmls(
  XMLS ${visit_readers_xmls})

# if you want to add vtk-module dependencies, simply add those in module.cmake,
# rather than explicitly adding those by using target_link_libraries.

# Required for paraDIS.
vtk_module_find_package(PRIVATE
  PACKAGE Boost)
vtk_module_link(ParaView::IOVisItBridge
  PRIVATE
    "$<BUILD_INTERFACE:Boost::boost>"
    ${visit_libraries})
