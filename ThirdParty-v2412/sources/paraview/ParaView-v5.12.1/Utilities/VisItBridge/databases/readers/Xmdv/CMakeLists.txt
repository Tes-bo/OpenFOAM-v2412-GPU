set(SOURCES
  avtXmdvFileFormat.C
  avtXmdvOptions.C
  avtXmdvWriter.C
  XmdvPluginInfo.C
  XmdvEnginePluginInfo.C
  XmdvCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItXmdvReader "1.0"
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtXmdvFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
