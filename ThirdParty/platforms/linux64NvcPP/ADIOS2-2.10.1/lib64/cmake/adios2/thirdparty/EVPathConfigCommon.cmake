include(CMakeFindDependencyMacro)

find_dependency(atl 2.2.1)
find_dependency(ffs 3.0.0)

if(NOT ON)
  find_dependency(dill 3.0.0)

  if()
    find_dependency(NVML)
  endif()

  if(NOT ON)
    if(TRUE)
      find_dependency(enet 1.3.14)
    endif()

    if(FALSE)
      find_dependency(LIBFABRIC)
    endif()

    if()
     find_dependency(NNTI)
    endif()
  endif()
endif()

include("${CMAKE_CURRENT_LIST_DIR}/EVPathConfigVersion.cmake")

include(FindPackageHandleStandardArgs)
set(EVPath_CONFIG "${CMAKE_CURRENT_LIST_FILE}")
find_package_handle_standard_args(EVPath CONFIG_MODE)

if(NOT TARGET EVPath::EVPath)
  include("${CMAKE_CURRENT_LIST_DIR}/EVPathTargets.cmake")
endif()

set(EVPath_LIBRARIES EVPath::EVPath)
set(EVPath_INCLUDE_DIRS
  $<TARGET_PROPERTY:EVPath::EVPath,INTERFACE_INCLUDE_DIRECTORIES>
)
