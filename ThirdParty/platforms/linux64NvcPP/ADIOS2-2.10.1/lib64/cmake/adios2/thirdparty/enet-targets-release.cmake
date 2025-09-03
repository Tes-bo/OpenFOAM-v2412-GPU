#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "enet::enet" for configuration "Release"
set_property(TARGET enet::enet APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(enet::enet PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libadios2_enet.so.2.10.1"
  IMPORTED_SONAME_RELEASE "libadios2_enet.so.2.10"
  )

list(APPEND _cmake_import_check_targets enet::enet )
list(APPEND _cmake_import_check_files_for_enet::enet "${_IMPORT_PREFIX}/lib64/libadios2_enet.so.2.10.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
