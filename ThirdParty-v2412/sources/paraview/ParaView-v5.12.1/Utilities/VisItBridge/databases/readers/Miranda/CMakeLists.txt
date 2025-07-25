set(SOURCES
  avtMirandaFileFormat.C
  avtMirandaOptions.C
  MirandaPluginInfo.C
  MirandaEnginePluginInfo.C
  MirandaCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItMirandaReader "1.0"
  VISIT_READER_TYPE "MTMD"
  VISIT_READER_NAME "avtMirandaFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
