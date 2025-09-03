#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "adios2::cxx11" for configuration "Release"
set_property(TARGET adios2::cxx11 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::cxx11 PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "adios2::core"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libadios2_cxx11.so.2.10.1"
  IMPORTED_SONAME_RELEASE "libadios2_cxx11.so.2.10"
  )

list(APPEND _cmake_import_check_targets adios2::cxx11 )
list(APPEND _cmake_import_check_files_for_adios2::cxx11 "${_IMPORT_PREFIX}/lib64/libadios2_cxx11.so.2.10.1" )

# Import target "adios2::cxx11_mpi" for configuration "Release"
set_property(TARGET adios2::cxx11_mpi APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(adios2::cxx11_mpi PROPERTIES
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "adios2::core_mpi"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib64/libadios2_cxx11_mpi.so.2.10.1"
  IMPORTED_SONAME_RELEASE "libadios2_cxx11_mpi.so.2.10"
  )

list(APPEND _cmake_import_check_targets adios2::cxx11_mpi )
list(APPEND _cmake_import_check_files_for_adios2::cxx11_mpi "${_IMPORT_PREFIX}/lib64/libadios2_cxx11_mpi.so.2.10.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
