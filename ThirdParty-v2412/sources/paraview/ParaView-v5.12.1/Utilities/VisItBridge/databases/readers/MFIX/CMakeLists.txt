set(SOURCES
  avtMFIXFileFormat.C
  avtMFIXOptions.C
  avtMFIXFileFormat.C
  MFIXPluginInfo.C
  MFIXEnginePluginInfo.C
  MFIXCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItMFIXReader "1.0"
  VISIT_READER_TYPE "MTMD"
  VISIT_READER_NAME "avtMFIXFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
