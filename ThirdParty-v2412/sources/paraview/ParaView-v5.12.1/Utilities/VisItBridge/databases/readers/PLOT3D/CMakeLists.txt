set(SOURCES
  avtPLOT3DFileFormat.C
  avtPLOT3DOptions.C
  PLOT3DCommonPluginInfo.C
  PLOT3DEnginePluginInfo.C
  PLOT3DFunctions.C
  PLOT3DPluginInfo.C
  vtkPLOT3DReader.C
  vtkPLOT3DReaderInternals.C
  )

ADD_VISIT_READER(VisItPLOT3DReader "1.0"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_NAME "avtPLOT3DFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
