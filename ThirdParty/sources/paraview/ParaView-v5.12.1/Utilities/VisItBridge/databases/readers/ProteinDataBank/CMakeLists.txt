SET(SOURCES
  avtProteinDataBankFileFormat.C
  avtProteinDataBankOptions.C
  ProteinDataBankPluginInfo.C
  ProteinDataBankEnginePluginInfo.C
  ProteinDataBankCommonPluginInfo.C
  )

ADD_VISIT_READER(VisItProteinDataBankReader "1.0"
  VISIT_READER_TYPE "STSD"
  VISIT_READER_NAME "avtProteinDataBankFileFormat"
  VISIT_READER_USES_OPTIONS ON
  SERVER_SOURCES ${SOURCES}
  )

