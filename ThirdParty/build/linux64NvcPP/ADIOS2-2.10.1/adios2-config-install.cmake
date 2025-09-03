cmake_minimum_required(VERSION 3.12)

set(_ADIOS2_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_LIST_DIR}")

if(NOT ON)
  if(NOT OFF)
    set(EVPath_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
  endif()

  if(NOT OFF)
    set(atl_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
  endif()

  if(NOT OFF)
    set(dill_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
  endif()

  if(NOT OFF)
    set(ffs_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
  endif()

  if(NOT OFF)
    set(enet_DIR ${CMAKE_CURRENT_LIST_DIR}/thirdparty)
  endif()
endif()

set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG "${CMAKE_CURRENT_LIST_FILE}")
include("${CMAKE_CURRENT_LIST_DIR}/adios2-config-common.cmake")

set(CMAKE_MODULE_PATH ${_ADIOS2_CMAKE_MODULE_PATH})
unset(_ADIOS2_CMAKE_MODULE_PATH)
