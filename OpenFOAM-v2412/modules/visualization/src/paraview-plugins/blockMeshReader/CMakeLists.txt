#------------------------------------------------------------------------------
project(ParaFoamBlockReader)

cmake_minimum_required(VERSION 3.8)
cmake_policy(SET CMP0002 NEW)   # For cmake >= 3
cmake_policy(SET CMP0057 NEW)   # For cmake >= 3.10.2
set(CMAKE_BUILD_TYPE Release)

# Fail if not building out-of-source
if (CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_CURRENT_BINARY_DIR)
    message(FATAL_ERROR
        "In-source builds disallowed. Use a separate build directory")
endif()

#### Not yet required
### MPI hint: prefer same MPI that OpenFOAM is using
## if (IS_DIRECTORY "$ENV{MPI_ARCH_PATH}")
##   set(ENV{MPI_HOME} "$ENV{MPI_ARCH_PATH}")
## endif()

# Take OSMESA hint from OpenFOAM MESA_ARCH_PATH
if (IS_DIRECTORY "$ENV{MESA_ARCH_PATH}")
  set(ENV{OSMESA_ROOT} "$ENV{MESA_ARCH_PATH}")
endif()

#-----------------------------------------------------------------------------
# Simple discovery and sanity checks

find_package(ParaView REQUIRED)

# Qualify OpenFOAM/ParaView libraries with -pvMAJ.MIN
# The pv_api can also be used for paraview-MAJ.MIN
if (DEFINED ENV{PARAVIEW_API})
    set(pv_api "$ENV{PARAVIEW_API}")
else()
    set(pv_api "${ParaView_VERSION_MAJOR}.${ParaView_VERSION_MINOR}")
    message("The PARAVIEW_API env was unset, using ${pv_api}")
endif()

# Pass through the OpenFOAM API value to the reader
set(foam_api "$ENV{FOAM_API}")
if (NOT "${foam_api}")
    set(foam_api "2.0")
endif()

set(foamPvLibQualifier "-pv${pv_api}")
set(PLUGIN_DIRNAME "paraview-${pv_api}")

#-----------------------------------------------------------------------------

message("================")
message("Using ParaView = $ENV{ParaView_DIR}")
message("OpenFOAM api   = ${foam_api}, libs = ${foamPvLibQualifier}")
## if (${MPI_C_FOUND})
## message("MPI include-dir= ${MPI_C_INCLUDE_DIRS}")
## endif()

#DEBUG get_cmake_property(_variableNames VARIABLES)
#DEBUG foreach (_variableName ${_variableNames})
#DEBUG     message(STATUS "${_variableName}=${${_variableName}}")
#DEBUG endforeach()

#-----------------------------------------------------------------------------
include(GNUInstallDirs)

# Default is PROJECT target, with supporting libraries from the same location
# Possible temporary override with 'STAGED_INSTALL_PREFIX' later

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    # Default install prefix same as FOAM_MODULE_PREFIX
    # (or from FOAM_USER_LIBBIN)
    if ("$ENV{FOAM_MODULE_PREFIX}")
        set(CMAKE_INSTALL_PREFIX "$ENV{FOAM_MODULE_PREFIX}"
            CACHE PATH "Install path prefix" FORCE)
    else()
        message(WARNING " FOAM_MODULE_PREFIX was unset, using FOAM_LIBBIN")
        get_filename_component(FOAM_MODULE_PREFIX $ENV{FOAM_LIBBIN} DIRECTORY)
        set(CMAKE_INSTALL_PREFIX "${FOAM_MODULE_PREFIX}"
            CACHE PATH "Install path prefix" FORCE)
    endif()
endif()
set(SUPPORTING_LIBBIN "${CMAKE_INSTALL_PREFIX}/lib")

message("Install prefix = ${CMAKE_INSTALL_PREFIX}")

include(CMakeLists-OpenFOAM.txt)
include(CMakeLists-library.txt)

#-----------------------------------------------------------------------------

# Configuration changes at paraview-5.7
set(OLD_PARAVIEW FALSE)
if (("${ParaView_VERSION_MAJOR}" LESS 5) OR
   (("${ParaView_VERSION_MAJOR}" EQUAL 5) AND ("${ParaView_VERSION_MINOR}" LESS 7)))
    set(OLD_PARAVIEW TRUE)
endif()

if ("${OLD_PARAVIEW}")
    message("================")
    include(${PARAVIEW_USE_FILE})
    include(CMakeLists-Project56.txt)
else()
    message("================")
    include(CMakeLists-Project.txt)
endif()

#-----------------------------------------------------------------------------
