include(../CMake/qtTestingMacroGenerateMocs.cmake)

IF(QtTesting_QT_VERSION VERSION_GREATER "5")
  FIND_PACKAGE(Qt6 REQUIRED QUIET COMPONENTS Core Test Widgets Gui)
  SET(TEST_LIBRARIES Qt6::Test Qt6::Core Qt6::Widgets Qt6::Gui)
ELSE()
  FIND_PACKAGE(Qt5 REQUIRED QUIET COMPONENTS Core Test Widgets Gui)
  SET(TEST_LIBRARIES Qt5::Test Qt5::Core Qt5::Widgets Qt5::Gui)
ENDIF()

set(KIT ${PROJECT_NAME})

set(TEST_SOURCES
  pqAbstractButtonEventTranslatorTest.cpp
  pqEventPlayerTest.cpp
  pqEventRecorderTest.cpp
  pqEventTranslatorTest.cpp
  pqDoubleSpinBoxEventPlayerTest.cpp
  pqDoubleSpinBoxEventTranslatorTest.cpp
  pqSpinBoxEventPlayerTest.cpp
  pqSpinBoxEventTranslatorTest.cpp
  pqTestUtilityTest.cpp
  )

set(TEST_MOC_HEADERS
  pqTest.h
  )

create_test_sourcelist(Tests ${KIT}CppTests.cxx
  ${TEST_SOURCES}
  )

set(TestsToRun ${Tests})
remove(TestsToRun ${KIT}CppTests.cxx)

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

if(QtTesting_QT_VERSION VERSION_GREATER "5")
  QT6_GENERATE_MOCS(${TEST_SOURCES})
  QT6_WRAP_CPP( TEST_MOC_SRCS ${TEST_MOC_HEADERS} )
else()
  QT5_GENERATE_MOCS(${TEST_SOURCES})
  QT5_WRAP_CPP( TEST_MOC_SRCS ${TEST_MOC_HEADERS} )
endif()

add_executable(${KIT}CppTests ${Tests} ${TEST_MOC_SRCS})
target_link_libraries(${KIT}CppTests ${PROJECT_NAME} ${TEST_LIBRARIES} qttesting)
set_target_properties(${KIT}CppTests PROPERTIES
  COMPILE_FLAGS "${Qt5Test_EXECUTABLE_COMPILE_FLAGS}")

macro(SIMPLE_TEST testname)
  add_test(NAME ${testname} COMMAND $<TARGET_FILE:${KIT}CppTests> ${testname} ${ARGN})
endmacro()

#
# Add Tests
#
SIMPLE_TEST( pqAbstractButtonEventTranslatorTest )
SIMPLE_TEST( pqEventPlayerTest )
SIMPLE_TEST( pqEventRecorderTest )
SIMPLE_TEST( pqDoubleSpinBoxEventPlayerTest )
SIMPLE_TEST( pqDoubleSpinBoxEventTranslatorTest )
SIMPLE_TEST( pqSpinBoxEventPlayerTest )
SIMPLE_TEST( pqSpinBoxEventTranslatorTest )
SIMPLE_TEST( pqEventTranslatorTest )
SIMPLE_TEST( pqTestUtilityTest )
