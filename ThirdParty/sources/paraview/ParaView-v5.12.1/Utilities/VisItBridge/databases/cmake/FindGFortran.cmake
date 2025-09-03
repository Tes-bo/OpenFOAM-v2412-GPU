set(_gfortran_extra_paths)
if (CMAKE_Fortran_COMPILER AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  get_filename_component(_gfortran_root "${CMAKE_Fortran_COMPILER}" DIRECTORY)
  get_filename_component(_gfortran_root "${_gfortran_root}" DIRECTORY)
  list(APPEND _gfortran_extra_paths
    "${_gfortran_root}/lib")
endif ()

find_library(GFortran_gfortran_LIBRARY
  NAMES gfortran
  HINTS
    ${_gfortran_extra_paths}
  PATH_SUFFIXES
    gcc/x86_64-linux-gnu/${CMAKE_Fortran_COMPILER_VERSION}/
    gcc/x86_64-redhat-linux/${CMAKE_Fortran_COMPILER_VERSION}/
  DOC "Path to the gfortran library")
mark_as_advanced(GFortran_gfortran_LIBRARY)

find_library(GFortran_quadmath_LIBRARY
  NAMES quadmath
  HINTS
    ${_gfortran_extra_paths}
  PATH_SUFFIXES
    gcc/x86_64-linux-gnu/${CMAKE_Fortran_COMPILER_VERSION}/
    gcc/x86_64-redhat-linux/${CMAKE_Fortran_COMPILER_VERSION}/
  DOC "Path to the quadmath library")
mark_as_advanced(GFortran_quadmath_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set GFortran_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GFortran DEFAULT_MSG GFortran_gfortran_LIBRARY GFortran_quadmath_LIBRARY)

if (GFortran_FOUND)
  set(GFortran_LIBRARIES
    "${GFortran_gfortran_LIBRARY}"
    "${GFortran_quadmatn_LIBRARY}")
  if (NOT TARGET GFortran::gfortran)
    add_library(GFortran::gfortran UNKNOWN IMPORTED)
    set_target_properties(GFortran::gfortran
      PROPERTIES
        IMPORTED_LOCATION "${GFortran_gfortran_LIBRARY}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran")
  endif ()
  if (NOT TARGET GFortran::quadmath)
    add_library(GFortran::quadmath UNKNOWN IMPORTED)
    set_target_properties(GFortran::quadmath
      PROPERTIES
        IMPORTED_LOCATION "${GFortran_quadmath_LIBRARY}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran")
  endif ()
  if (NOT TARGET GFortran::GFortran)
    add_library(GFortran::GFortran INTERFACE IMPORTED)
    set_target_properties(GFortran::GFortran
      PROPERTIES
        INTERFACE_LINK_LIBRARIES "GFortran::gfortran;GFortran::quadmath")
  endif ()
endif ()
