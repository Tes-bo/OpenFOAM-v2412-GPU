find_path(Mili_INCLUDE_DIR
  NAMES mili.h
  PATH_SUFFIXES mili
  DOC   "Path to the Mili headers")
mark_as_advanced(Mili_INCLUDE_DIR)

find_library(Mili_LIBRARY
  NAMES mili
  DOC   "Path to the Mili library")
mark_as_advanced(Mili_LIBRARY)

# handle the QUIETLY and REQUIRED arguments and set Mili_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Mili DEFAULT_MSG Mili_LIBRARY Mili_INCLUDE_DIR)

if (Mili_FOUND)
  set(Mili_LIBRARIES "${Mili_LIBRARY}")
  set(Mili_INCLUDE_DIRS "${Mili_INCLUDE_DIR}")
  if (NOT TARGET Mili::Mili)
    add_library(Mili::Mili UNKNOWN IMPORTED)
    set_target_properties(Mili::Mili
      PROPERTIES
        IMPORTED_LOCATION "${Mili_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${Mili_INCLUDE_DIR}")
  endif ()
endif ()
