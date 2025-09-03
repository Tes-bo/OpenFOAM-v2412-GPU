set(SOURCES
  avtDyna3DFileFormat.C
  Dyna3DFile.cpp
  MaterialProperties.cpp
  )

ADD_VISIT_READER(VisItDyna3DReader "1.0"
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtDyna3DFileFormat"
  SERVER_SOURCES ${SOURCES}
  )
