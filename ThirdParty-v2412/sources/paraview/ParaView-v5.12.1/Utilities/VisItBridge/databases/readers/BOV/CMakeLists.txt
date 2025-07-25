set(SOURCES
  avtBOVFileFormat.C
  avtBOVOptions.C
  avtBOVWriter.C
  BOVPluginInfo.C
  BOVEnginePluginInfo.C
  BOVCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItBOVReader "1.0"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_NAME "avtBOVFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
