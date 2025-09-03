find_path(SILO_INCLUDE_DIR
  NAMES silo.h
  DOC   "Path to the SILO include directory")
mark_as_advanced(SILO_INCLUDE_DIR)

find_library(SILO_LIBRARY
  NAMES silo siloh5 silohdf5
  DOC   "Path to the SILO library")
mark_as_advanced(SILO_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set SILO_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SILO DEFAULT_MSG SILO_LIBRARY SILO_INCLUDE_DIR)

if (SILO_FOUND)
  set(SILO_LIBRARIES "${SILO_LIBRARY}")
  set(SILO_INCLUDE_DIRS "${SILO_INCLUDE_DIR}")
  if (NOT TARGET SILO::SILO)
    add_library(SILO::SILO UNKNOWN IMPORTED)
    set_target_properties(SILO::SILO
      PROPERTIES
        IMPORTED_LOCATION "${SILO_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${SILO_INCLUDE_DIR}")
  endif ()
endif ()
