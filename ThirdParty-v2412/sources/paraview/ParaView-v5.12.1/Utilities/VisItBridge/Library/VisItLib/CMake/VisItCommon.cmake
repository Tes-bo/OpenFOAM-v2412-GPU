#Setup the version since some files need this
set(VISIT_VERSION "2.7.0")

if (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
endif (COMMAND cmake_policy)

#disable compiler warnings from the bridge
if (FALSE)
option(VISIT_DISABLE_COMPILER_WARNINGS "Disable compiler warnings" ON)
mark_as_advanced(VISIT_DISABLE_COMPILER_WARNINGS)
else ()
set(VISIT_DISABLE_COMPILER_WARNINGS ON)
endif ()
if(VISIT_DISABLE_COMPILER_WARNINGS)
  if(WIN32)
    if (MSVC)
      string(REGEX REPLACE "(^| )([/-])W[0-9]( |$)" " "
          CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
      set(CMAKE_C_FLAGS "${CMAKE_CXX_FLAGS} /w")

      string(REGEX REPLACE "(^| )([/-])W[0-9]( |$)" " "
          CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /w")
    endif(MSVC)
  else(WIN32)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -w")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
  endif(WIN32)
endif(VISIT_DISABLE_COMPILER_WARNINGS)

#set up some vars we need to compile
set(VISIT_STATIC)
if (BUILD_SHARED_LIBS)
  set(VISIT_STATIC 0)
  add_definitions(-DVISIT_BUILD_SHARED_LIBS)
else(VISIT_STATIC)
  set(VISIT_STATIC 1)
  add_definitions(-DVISIT_STATIC -DGLEW_STATIC)
endif()

set(VISIT_SOURCE_DIR ${VisItBridgePlugin_SOURCE_DIR})
set(VISIT_BINARY_DIR ${VisItBridgePlugin_BINARY_DIR})
set(VISIT_CMAKE_DIR ${VISIT_SOURCE_DIR}/CMake )

#include the visit cmake directory on the cmake search path
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${VISIT_CMAKE_DIR})

#include the visit install target and plugin function
MACRO(VISIT_INSTALL_TARGETS target)
  _vtk_module_apply_properties(${target})
  _vtk_module_install(${target})
ENDMACRO(VISIT_INSTALL_TARGETS)

FUNCTION(ADD_PARALLEL_LIBRARY target)
    ADD_LIBRARY(${target} ${ARGN})
    IF(VISIT_PARALLEL_CXXFLAGS)
        SET(PAR_COMPILE_FLAGS "")
        FOREACH (X ${VISIT_PARALLEL_CXXFLAGS})
            SET(PAR_COMPILE_FLAGS "${PAR_COMPILE_FLAGS} ${X}")
        ENDFOREACH(X)
        SET_TARGET_PROPERTIES(${target} PROPERTIES
            COMPILE_FLAGS ${PAR_COMPILE_FLAGS}
        )
        IF(VISIT_PARALLEL_LINKER_FLAGS)
            SET(PAR_LINK_FLAGS "")
            FOREACH(X ${VISIT_PARALLEL_LINKER_FLAGS})
                SET(PAR_LINK_FLAGS "${PAR_LINK_FLAGS} ${X}")
            ENDFOREACH(X)
            SET_TARGET_PROPERTIES(${target} PROPERTIES
                LINK_FLAGS ${PAR_LINK_FLAGS}
            )
        ENDIF(VISIT_PARALLEL_LINKER_FLAGS)

        IF(VISIT_PARALLEL_RPATH)
            SET(PAR_RPATHS "")
            FOREACH(X ${CMAKE_INSTALL_RPATH})
                SET(PAR_RPATHS "${PAR_RPATHS} ${X}")
            ENDFOREACH(X)
            FOREACH(X ${VISIT_PARALLEL_RPATH})
                SET(PAR_RPATHS "${PAR_RPATHS} ${X}")
            ENDFOREACH(X)
            SET_TARGET_PROPERTIES(${target} PROPERTIES
                INSTALL_RPATH ${PAR_RPATHS}
            )
        ENDIF(VISIT_PARALLEL_RPATH)

        IF(NOT VISIT_NOLINK_MPI_WITH_LIBRARIES)
            target_link_libraries(${target} ${VISIT_PARALLEL_LIBS})
        ENDIF(NOT VISIT_NOLINK_MPI_WITH_LIBRARIES)
    ENDIF(VISIT_PARALLEL_CXXFLAGS)
ENDFUNCTION(ADD_PARALLEL_LIBRARY)

#set up MPI
if (TARGET VTK::mpi)
  set(VISIT_PARALLEL 1)
  set(VISIT_PARALLEL_LIBS VTK::mpi)
else ()
  set(VISIT_PARALLEL 0)
  set(VISIT_PARALLEL_LIBS VTK::mpi)
endif ()

# setup to use vtkzlib
set(ZLIB_LIB VTK::zlib)

#block out most of the warnings in visit on windows
if (WIN32)
  add_definitions(-D_USE_MATH_DEFINES)
  if (MSVC_VERSION GREATER 1400)
      add_definitions(-D_CRT_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_DEPRECATE)
      add_definitions(-D_SCL_NO_DEPRECATE -D_SCL_SECURE_NO_DEPRECATE)
      add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  endif (MSVC_VERSION GREATER 1400)
endif(WIN32)

#-----------------------------------------------------------------------------
# Setup ParaView and Common includes before the visit-config.h so that
# we can use utilities like HDF5 from ParaView
#-----------------------------------------------------------------------------


# Set up easy to use includes for the common directory
set(VISIT_COMMON_INCLUDES
    ${VISIT_BINARY_DIR}
    ${VISIT_SOURCE_DIR}/common/Exceptions/Database
    ${VISIT_SOURCE_DIR}/common/Exceptions/Pipeline
    ${VISIT_SOURCE_DIR}/common/Exceptions/Plotter
    ${VISIT_SOURCE_DIR}/common/comm
    ${VISIT_SOURCE_DIR}/common/expr
    ${VISIT_SOURCE_DIR}/common/misc
    ${VISIT_SOURCE_DIR}/common/parser
    ${VISIT_SOURCE_DIR}/common/plugin
    ${VISIT_SOURCE_DIR}/common/proxybase
    ${VISIT_SOURCE_DIR}/common/state
    ${VISIT_SOURCE_DIR}/common/utility
  )

# Watch out, you need to make sure common/parser is always in front of
# python2.X includes
include_directories(BEFORE ${VISIT_COMMON_INCLUDES})

#-----------------------------------------------------------------------------
# Setup Vars for visit-config.h
#-----------------------------------------------------------------------------


set(VISIT_DBIO_ONLY ON) #Build only visitconvert and engine plugins
if(VISIT_DBIO_ONLY)
    add_definitions(-DDBIO_ONLY)
endif(VISIT_DBIO_ONLY)

#Check to see if ofstreams rdbuf is public. If it is NOT public set NO_SETBUF
try_compile(VISIT_COMPILER_FSTREAM_WORKAROUND
        ${CMAKE_CURRENT_BINARY_DIR}/CMakeTmp
        ${VISIT_CMAKE_DIR}/testFStream.cxx)
if (NOT VISIT_COMPILER_FSTREAM_WORKAROUND)
   set(NO_SETBUF 1)
endif (NOT VISIT_COMPILER_FSTREAM_WORKAROUND)

# Set the slash characters based on the platform
if(WIN32)
    set(VISIT_SLASH_CHAR   "'\\\\'")
    set(VISIT_SLASH_STRING "\"\\\\\"")
else(WIN32)
    set(VISIT_SLASH_CHAR   "'/'")
    set(VISIT_SLASH_STRING "\"/\"")
endif(WIN32)

# Check for plugin extension
if(VISIT_STATIC)
  if(WIN32)
    set(VISIT_PLUGIN_EXTENSION   ".dll")
  else(WIN32)
    set(VISIT_PLUGIN_EXTENSION   ".lib")
  endif(WIN32)
else(VISIT_STATIC)
    if(WIN32)
        set(VISIT_PLUGIN_EXTENSION   ".dll")
    else(WIN32)
        if(APPLE)
            set(VISIT_PLUGIN_EXTENSION   ".dylib")
        else(APPLE)
            set(VISIT_PLUGIN_EXTENSION   ".so")
        endif(APPLE)
    endif(WIN32)
endif(VISIT_STATIC)

#-----------------------------------------------------------------------------
# Setup lib settings
#-----------------------------------------------------------------------------

#setup non third party vtk utilities
set(HAVE_LIBHDF5 1)
set(HAVE_LIBNETCDF 1)
set(HAVE_NETCDF_H 1)
set(HAVE_LIBEXODUSII 1)

#-----------------------------------------------------------------------------
# Detect packages here. We could probably write macros that we can include from
# elsewhere for this.
#-----------------------------------------------------------------------------
include(CheckIncludeFiles)
include(CMakeBackwardCompatibilityC)
include(CMakeBackwardCompatibilityCXX)
include(CheckTypeSize)
include(CheckFunctionExists)
include(CheckSymbolExists)
include(TestBigEndian)
#include(FindOpenGL)

check_include_files (fcntl.h     HAVE_FCNTL_H)
check_include_files (inttypes.h  HAVE_INTTYPES_H)
check_include_files (malloc.h    HAVE_MALLOC_H)
check_include_files (limits.h    HAVE_LIMITS_H)
check_include_files (memory.h    HAVE_MEMORY_H)
check_include_files (stdint.h    HAVE_STDINT_H)
check_include_files (stdlib.h    HAVE_STDLIB_H)
check_include_files (strings.h   HAVE_STRINGS_H)
check_include_files (string.h    HAVE_STRING_H)
check_include_files (sys/time.h  HAVE_SYS_TIME_H)
check_include_files (sys/types.h HAVE_SYS_TYPES_H)
check_include_files (sys/stat.h  HAVE_SYS_STAT_H)
check_include_files (unistd.h    HAVE_UNISTD_H)
check_include_files (stdbool.h   HAVE_STDBOOL_H)

# Check for type sizes, endian
set(SIZEOF_BOOLEAN              ${CMAKE_SIZEOF_BOOLEAN})
set(SIZEOF_CHAR                 ${CMAKE_SIZEOF_CHAR})
set(SIZEOF_DOUBLE               ${CMAKE_SIZEOF_DOUBLE})
set(SIZEOF_FLOAT                ${CMAKE_SIZEOF_FLOAT})
set(SIZEOF_INT                  ${CMAKE_SIZEOF_INT})
set(SIZEOF_LONG                 ${CMAKE_SIZEOF_LONG})
set(SIZEOF_LONG_DOUBLE          ${CMAKE_SIZEOF_LONG_DOUBLE})
set(SIZEOF_LONG_FLOAT           ${CMAKE_SIZEOF_LONG_FLOAT})
set(SIZEOF_LONG_LONG            ${CMAKE_SIZEOF_LONG_LONG})
set(SIZEOF_SHORT                ${CMAKE_SIZEOF_SHORT})
set(SIZEOF_UNSIGNED_CHAR        ${CMAKE_SIZEOF_UNSIGNED_CHAR})
set(SIZEOF_UNSIGNED_INT         ${CMAKE_SIZEOF_UNSIGNED_INT})
set(SIZEOF_UNSIGNED_LONG        ${CMAKE_SIZEOF_UNSIGNED_LONG})
set(SIZEOF_UNSIGNED_LONG_LONG   ${CMAKE_SIZEOF_UNSIGNED_LONG_LONG})
set(SIZEOF_UNSIGNED_SHORT       ${CMAKE_SIZEOF_UNSIGNED_SHORT})
set(SIZEOF_VOID_P               ${CMAKE_SIZEOF_VOID_P})
set(CMAKE_REQUIRED_DEFINITIONS -D_LARGEFILE64_SOURCE)
check_type_size("off64_t" SIZEOF_OFF64_T)
test_big_endian(WORDS_BIGENDIAN)

# manually check for socklen_t as CHECK_SYMBOL_EXISTS
# doesn't appear to work on linux (at least)
try_compile(HAVE_SOCKLEN_T
    ${CMAKE_CURRENT_BINARY_DIR}
    ${VISIT_SOURCE_DIR}/CMake/TestSocklenT.c
    OUTPUT_VARIABLE SLT
)
if (HAVE_SOCKLEN_T)
    set(HAVE_SOCKLEN_T 1 CACHE INTERNAL "support for socklen_t")
else(HAVE_SOCKLEN_T)
    set(HAVE_SOCKLEN_T 0 CACHE INTERNAL "support for socklen_t")
endif (HAVE_SOCKLEN_T)
