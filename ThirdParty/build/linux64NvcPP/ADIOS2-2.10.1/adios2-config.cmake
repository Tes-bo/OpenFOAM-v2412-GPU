set(_ADIOS2_CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})
list(INSERT CMAKE_MODULE_PATH 0 "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake")

if(NOT ON)
  set(atl_DIR /home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/thirdparty/atl/atl)
  set(dill_DIR /home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/thirdparty/dill/dill)
  set(ffs_DIR /home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/thirdparty/ffs/ffs)
endif()

if(TRUE)
  set(EVPath_DIR /home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/thirdparty/EVPath/EVPath)
  if(NOT ON)
    set(enet_DIR /home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/thirdparty/enet/enet)
  endif()
endif()

set(${CMAKE_FIND_PACKAGE_NAME}_CONFIG "${CMAKE_CURRENT_LIST_FILE}")
include("${CMAKE_CURRENT_LIST_DIR}/adios2-config-common.cmake")

set(CMAKE_MODULE_PATH ${_ADIOS2_CMAKE_MODULE_PATH})
unset(_ADIOS2_CMAKE_MODULE_PATH)
