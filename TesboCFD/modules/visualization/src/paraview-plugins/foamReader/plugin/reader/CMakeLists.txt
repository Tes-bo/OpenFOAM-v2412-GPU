#-----------------------------------------------------------------------------
set(classes vtkPVFoamReader)

set(TARGET_NAME FoamReader)

vtk_module_add_module(${TARGET_NAME}
  FORCE_STATIC
  CLASSES ${classes})

# Predefined values from CMakeLists-library.txt
target_include_directories(${TARGET_NAME}
    PUBLIC
    ${OPENFOAM_INCLUDE_DIRS}
)


# With cmake 3.12 and newer, can simply use target_link_directories
# ----
# target_link_directories(${TARGET_NAME}
#     PUBLIC
#     ${FOAM_LIBBIN}
#     ${SUPPORTING_LIBBIN}
# )

# target_link_libraries(${TARGET_NAME}
#     PUBLIC
#     ${OPENFOAM_LIBNAMES}
#     ${SUPPORTING_LIBNAMES}
# )
# ----


# Older cmake. Do manually
# ----

# Predefined values from CMakeLists-library.txt
link_directories(${OPENFOAM_LINK_DIRS})

unset(OPENFOAM_LIBS)
foreach(libname ${OPENFOAM_LIBNAMES})
    unset(libpath CACHE)
    find_library(libpath ${libname} PATHS ${FOAM_LIBBIN})
    list(APPEND OPENFOAM_LIBS ${libpath})
endforeach(libname)

unset(SUPPORTING_LIBS)
foreach(libname ${SUPPORTING_LIBNAMES})
    unset(libpath CACHE)
    find_library(libpath ${libname} PATHS ${SUPPORTING_LIBBIN})
    list(APPEND SUPPORTING_LIBS ${libpath})
endforeach(libname)

## message("Uses libs ${SUPPORTING_LIBNAMES}")

target_link_libraries(${TARGET_NAME}
    PUBLIC
    Qt5::Core
    Qt5::Gui
    ${OPENFOAM_LIBS}
    ${SUPPORTING_LIBS}
)

#-----------------------------------------------------------------------------
