SET(SOURCES
  avtChomboFileFormat.C
  avtChomboOptions.C
  ChomboPluginInfo.C
  ChomboEnginePluginInfo.C
  ChomboCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItChomboReader "1.0"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_NAME "avtChomboFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
