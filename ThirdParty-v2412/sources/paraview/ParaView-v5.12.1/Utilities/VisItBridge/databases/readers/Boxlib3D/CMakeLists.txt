SET(SOURCES
  avtBoxlib3DFileFormat.C
  BoxLibHelper.C
  )

if(Boxlib_FOUND)
  ADD_VISIT_READER(VisItBoxlib3DReader "1.0"
    VISIT_READER_TYPE "STMD"
    VISIT_READER_NAME "avtBoxlib3DFileFormat"
    VISIT_READER_USES_OPTIONS OFF
    SERVER_SOURCES ${SOURCES}
    )
endif()
