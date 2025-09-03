# Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
# Project developers.  See the top-level LICENSE file for dates and other
# details.  No copyright assignment is required to contribute to VisIt.

#*****************************************************************************
# Modifications:
#
#****************************************************************************

SET(AVTMATH_SOURCES
#avtFactorial.C
#avtHermitePolynomial.C
avtMatrix.C
#avtPlane.C
avtPolynomial.C
#avtQuaternion.C
#avtTrackball.C
#avtVecArray.C
avtVector.C
)

if (FALSE)
INCLUDE_DIRECTORIES(
${CMAKE_CURRENT_SOURCE_DIR}
${VISIT_COMMON_INCLUDES}
${VTK_INCLUDE_DIRS}
)

# Add link directories
LINK_DIRECTORIES(${LIBRARY_OUTPUT_DIRECTORY} )
endif ()

ADD_LIBRARY(avtmath ${AVTMATH_SOURCES})
#TARGET_LINK_LIBRARIES(avtmath visitcommon vtkCommonCore)
TARGET_LINK_LIBRARIES(avtmath PRIVATE visitcommon VTK::CommonCore)
target_include_directories(avtmath
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

VISIT_INSTALL_TARGETS(avtmath)
