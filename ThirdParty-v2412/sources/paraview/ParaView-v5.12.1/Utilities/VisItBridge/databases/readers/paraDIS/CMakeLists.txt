set(SOURCES
  avtparaDISFileFormat.C
  avtparaDISOptions.C
  paraDISCommonPluginInfo.C
  paraDISEnginePluginInfo.C
  paraDISPluginInfo.C
  Dumpfile.C
  parallelParaDIS.C
  paraDIS_lib/paradis.C
  paraDIS_lib/paradis_c_interface.C
  RC_cpp_lib/RC_c_lib/args.c
  RC_cpp_lib/RC_c_lib/debugutil.c
  RC_cpp_lib/RC_c_lib/fileutils.c
  RC_cpp_lib/RC_c_lib/inventor.c
  RC_cpp_lib/RangeList.C
  RC_cpp_lib/timer.C
  RC_cpp_lib/pathutil.C
  )

ADD_VISIT_READER(VisItParaDISReader "1.0"
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtparaDISFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )
