#------------------------------------------------------------------------------
project(runTimePostProcessing)

cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0002 NEW)   # For cmake >= 3
cmake_policy(SET CMP0057 NEW)   # For cmake >= 3.10.2
set(CMAKE_BUILD_TYPE Release)

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

# MPI hint: prefer same MPI that OpenFOAM is using
if (IS_DIRECTORY "$ENV{MPI_ARCH_PATH}")
  set(ENV{MPI_HOME} "$ENV{MPI_ARCH_PATH}")
endif()

# Take OSMESA hint from OpenFOAM MESA_ARCH_PATH
if (IS_DIRECTORY "$ENV{MESA_ARCH_PATH}")
  set(ENV{OSMESA_ROOT} "$ENV{MESA_ARCH_PATH}")
endif()

#-----------------------------------------------------------------------------
# Simple discovery and sanity checks

set(ParaView_VERSION 0)
set(VTK_VERSION 0)

if (EXISTS "$ENV{VTK_DIR}")
    message("Building with VTK from $ENV{VTK_DIR}")
    find_package(VTK REQUIRED)
elseif (EXISTS "$ENV{ParaView_DIR}")
    message("Building with Paraview from $ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED)
else()
    message(FATAL_ERROR "VTK not found using VTK_DIR or ParaView_DIR")
endif()

#-----------------------------------------------------------------------------
# Earlier versions without modules and/or VTK_VERSION
#
# ParaView = 5.6.3 (VTK equiv 8.2.0)
# ParaView = 5.7.0 (VTK 8.90.0)
# ParaView = 5.8.1 (VTK 8.90.0)
# ParaView = 5.9.0 (VTK 9.0.1)

# VTK has modules yes/no.
# Added in VTK9 but also in VTK-8.90 (paraview flavour for 5.7, 5.8)
set(OLDVTK_NO_MODULES FALSE)
if (("${VTK_MAJOR_VERSION}" LESS 8) OR
   (("${VTK_MAJOR_VERSION}" EQUAL 8) AND ("${VTK_MINOR_VERSION}" LESS 10)))
    set(OLDVTK_NO_MODULES TRUE)
endif()

if (VTK_MAJOR_VERSION AND NOT VTK_VERSION)
  set(VTK_VERSION
      "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}.${VTK_BUILD_VERSION}")
  message("Synthesized VTK version: ${VTK_VERSION}")
endif()

#-----------------------------------------------------------------------------

# New ParaView/VTK: load optional MPI via VTK
# Old ParaView/VTK: has VTK_USE_FILE, manually detect ParallelMPI

if ("${OLDVTK_NO_MODULES}")
    # No modules etc - use old method for loading things
    if (VTK_USE_FILE)
        message("VTK from ${VTK_USE_FILE}")
        include(${VTK_USE_FILE})
    endif()
else()
   # Load optional components via modules
   if (${ParaView_VERSION})
       find_package(ParaView REQUIRED
           OPTIONAL_COMPONENTS VTK::ParallelMPI)
   else()
       find_package(VTK REQUIRED
           OPTIONAL_COMPONENTS ParallelMPI)
   endif()
endif()

#-----------------------------------------------------------------------------
include(GNUInstallDirs)

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
if (${ParaView_VERSION})
message("Using ParaView = ${ParaView_VERSION} (VTK ${VTK_VERSION})")
else()
message("Using VTK ${VTK_VERSION}")
endif()
if (${MPI_C_FOUND})
message("MPI include-dir= ${MPI_C_INCLUDE_DIRS}")
endif()

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-library.txt)

# Default is PROJECT target, with supporting libraries from FOAM_LIBBIN
# Possible temporary override with 'STAGED_INSTALL_PREFIX' later
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    get_filename_component(FOAM_PLATFORMS ${FOAM_LIBBIN} DIRECTORY)
    set(CMAKE_INSTALL_PREFIX "${FOAM_PLATFORMS}"
        CACHE PATH "Install path prefix" FORCE)
endif()
message("Install prefix = ${CMAKE_INSTALL_PREFIX}")


message("================")
include(CMakeLists-Project.txt)

#-----------------------------------------------------------------------------
