SET(SOURCES
  avtCurve2DFileFormat.C
  avtCurve2DOptions.C
  avtCurve2DWriter.C
  Curve2DPluginInfo.C
  Curve2DEnginePluginInfo.C
  Curve2DCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItCurve2DReader "1.0"
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtCurve2DFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
