# Install script for directory: /home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/install/pre

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/sibo/TesboCFD/ThirdParty-v2412/platforms/linux64NvcPP/ADIOS2-2.10.1")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "adios2_core-development" OR NOT CMAKE_INSTALL_COMPONENT)
  
  message("Pre-installation cleanup of CMake files")
  file(REMOVE_RECURSE "/home/sibo/TesboCFD/ThirdParty-v2412/platforms/linux64NvcPP/ADIOS2-2.10.1/lib64/cmake/adios2")

endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "adios2_core-development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/adios2" TYPE FILE FILES
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/CMakeFindDependencyMacro.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindBZip2.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindCrayDRC.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindDAOS.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindDataSpaces.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindHDF5.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindIME.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindLIBFABRIC.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindMPI.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindPkgConfig.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindPython.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindPythonModule.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindSZ.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindSodium.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindUCX.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/FindZeroMQ.cmake"
    "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/Findpugixml.cmake"
    )
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "adios2_core-development" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib64/cmake/adios2/upstream" TYPE DIRECTORY FILES "/home/sibo/TesboCFD/ThirdParty-v2412/sources/adios/ADIOS2-2.10.1/cmake/upstream/")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/home/sibo/TesboCFD/ThirdParty-v2412/build/linux64NvcPP/ADIOS2-2.10.1/cmake/install/pre/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
