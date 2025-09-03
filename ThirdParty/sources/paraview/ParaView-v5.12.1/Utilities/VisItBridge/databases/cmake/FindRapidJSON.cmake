set(RapidJSON_FOUND 0)

set(_RapidJSON_quiet_arg)
if (RapidJSON_FIND_QUIETLY)
  set(_RapidJSON_quiet_arg QUIET)
endif ()
find_package(RapidJSON CONFIG ${_RapidJSON_quiet_arg})
unset(_RapidJSON_quiet_arg)

if (NOT RapidJSON_FOUND)
  find_file(RapidJSON_INCLUDE_DIR
    NAMES rapidjson/rapidjson.h
    DOC "RapidJSON include directory")
  mark_as_advanced(RapidJSON_INCLUDE_DIR)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(RapidJson
    REQUIRED_VARS RapidJSON_INCLUDE_DIR)

  if (RapidJSON_FOUND)
    set(RapidJSON_INCLUDE_DIRS "${RapidJSON_INCLUDE_DIR}")
  endif ()
endif ()

# Upstream added the target in March 2020.
if (RapidJSON_FOUND)
  if (NOT TARGET rapidjson)
    add_library(rapidjson INTERFACE IMPORTED)
    set_property(TARGET rapidjson
      PROPERTY
        INTERFACE_INCLUDE_DIRECTORIES "${RapidJSON_INCLUDE_DIRS}")
  endif ()
endif ()
