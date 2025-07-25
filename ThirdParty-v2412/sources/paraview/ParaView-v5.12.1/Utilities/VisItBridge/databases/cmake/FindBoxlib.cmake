find_path(Boxlib_INCLUDE_DIR
  NAMES BoxLib.H
  DOC   "Path to the Boxlib headers")
mark_as_advanced(Boxlib_INCLUDE_DIR)

find_library(Boxlib_C_LIBRARY
  NAMES cboxlib
  DOC   "Path to the Boxlib C library")
mark_as_advanced(Boxlib_C_LIBRARY)

find_library(Boxlib_Fortran_LIBRARY
  NAMES fboxlib
  DOC   "Path to the Boxlib Fortran library")
mark_as_advanced(Boxlib_Fortran_LIBRARY)

if (Boxlib_FOUND)
  find_package(MPI QUIET COMPONENTS C Fortran)
  if (NOT MPI_FOUND)
    set(Boxlib_FOUND OFF)
    string(APPEND Boxlib_NOT_FOUND_MESSAGE
      "Boxlib needs C and Fortran MPI bindings.")
  endif ()
endif ()

if (Boxlib_FOUND)
  find_package(GFortran QUIET)
  if (NOT GFortran_FOUND)
    set(Boxlib_FOUND OFF)
    string(APPEND Boxlib_NOT_FOUND_MESSAGE
      "Boxlib needs GFortran.")
  endif ()
endif ()

# handle the QUIETLY and REQUIRED arguments and set Boxlib_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Boxlib DEFAULT_MSG Boxlib_C_LIBRARY Boxlib_Fortran_LIBRARY Boxlib_INCLUDE_DIR)

if (Boxlib_FOUND)
  set(Boxlib_LIBRARIES
    "${Boxlib_C_LIBRARY}"
    "${Boxlib_Fortran_LIBRARY}")
  set(Boxlib_INCLUDE_DIRS "${Boxlib_INCLUDE_DIR}")
  if (NOT TARGET Boxlib::C)
    add_library(Boxlib::C UNKNOWN IMPORTED)
    set_target_properties(Boxlib::C
      PROPERTIES
        IMPORTED_LOCATION "${Boxlib_C_LIBRARY}"
        INTERFACE_LINK_LIBRARIES "MPI::C"
        INTERFACE_INCLUDE_DIRECTORIES "${Boxlib_INCLUDE_DIR}")
  endif ()
  if (NOT TARGET Boxlib::Fortran)
    add_library(Boxlib::Fortran UNKNOWN IMPORTED)
    set_target_properties(Boxlib::Fortran
      PROPERTIES
        IMPORTED_LOCATION "${Boxlib_Fortran_LIBRARY}"
        IMPORTED_LINK_INTERFACE_LANGUAGES "Fortran"
        INTERFACE_LINK_LIBRARIES "MPI::Fortran;GFortran::GFortran"
        INTERFACE_INCLUDE_DIRECTORIES "${Boxlib_INCLUDE_DIR}")
  endif ()
  if (NOT TARGET Boxlib::Boxlib)
    add_library(Boxlib::Boxlib INTERFACE IMPORTED)
    set_target_properties(Boxlib::Boxlib
      PROPERTIES
        INTERFACE_LINK_LIBRARIES "Boxlib::C;Boxlib::Fortran")
  endif ()
endif ()
