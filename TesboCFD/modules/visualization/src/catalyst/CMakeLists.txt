#------------------------------------------------------------------------------
project(catalyst)

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

#------------------------------------------------------------------------------
# Simple discovery and sanity checks

if (EXISTS "$ENV{ParaView_DIR}")
    find_package(ParaView REQUIRED)
else()
    message(FATAL_ERROR "ParaView_DIR not set")
endif()

# Define pv_api (MAJ.MIN), and pv_version (MAJ.MIN.PATCH)
set(pv_api "${ParaView_VERSION_MAJOR}.${ParaView_VERSION_MINOR}")
set(pv_version "${ParaView_VERSION}")

# VTK has modules yes/no.
# Added in VTK9 but also in VTK-8.90 (paraview flavour for 5.7, 5.8)
set(OLDVTK_NO_MODULES FALSE)
if (("${VTK_MAJOR_VERSION}" LESS 8) OR
   (("${VTK_MAJOR_VERSION}" EQUAL 8) AND ("${VTK_MINOR_VERSION}" LESS 10)))
    set(OLDVTK_NO_MODULES TRUE)
    message("Using older VTK (no modules)")
endif()

#-----------------------------------------------------------------------------

# New ParaView/VTK: load optional MPI via VTK
# Old ParaView/VTK: has VTK_USE_FILE, manually detect ParallelMPI

unset(PROJECT_CATALYST_COMPONENTS)

if ("${OLDVTK_NO_MODULES}")
    # No modules etc - use old method for loading things
    find_package(ParaView REQUIRED
        COMPONENTS vtkPVPythonCatalyst
        OPTIONAL_COMPONENTS vtkParallelMPI)
    include(${PARAVIEW_USE_FILE})
else()
    set(PROJECT_CATALYST_COMPONENTS
        ParaView::PythonCatalyst
        ParaView::RemotingCore
    )

    find_package(ParaView REQUIRED
        COMPONENTS PythonCatalyst RemotingCore
        OPTIONAL_COMPONENTS VTK::ParallelMPI)

    # Catalyst package (ParaView > 5.8)
    if (("${PARAVIEW_VERSION_MAJOR}" GREATER 5) OR
        ("${PARAVIEW_VERSION_MINOR}" GREATER 8))
        find_package(catalyst REQUIRED)
    endif()
endif()

#-----------------------------------------------------------------------------
include(GNUInstallDirs)

message("================")
message("Using OpenFOAM = $ENV{WM_PROJECT_DIR}")
message("Using ParaView = $ENV{ParaView_DIR}")
if (${MPI_C_FOUND})
message("MPI include-dir= ${MPI_C_INCLUDE_DIRS}")
endif()

include(CMakeLists-OpenFOAM.txt)

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
