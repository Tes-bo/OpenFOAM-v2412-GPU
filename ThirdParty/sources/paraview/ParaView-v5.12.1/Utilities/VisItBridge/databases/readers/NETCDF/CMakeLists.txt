PROJECT(VisItNETCDF)

#############
## WARNING ##
#############
#VisIt interface framework is really weird. This is a quick first
#pass at figuring out a semi easy way to emulate it in ParaView world

#macro has to return true/false ( users needs to write it)
#this macro is called from the CanReadFile method of each class
set(INTERFACE_MACRO_CALL "ValidNETCDFFile(fname)")
set(INTERFACE_MACRO_FILE "ValidNETCDFFile") #needs to be a .h

SET(SOURCES
  avtBasicNETCDFFileFormat.C
  avtADAPTFileFormat.C
  avtLODIFileFormat.C
  avtLODIParticleFileFormat.C
  avtFVCOM_STSDFileFormat.C
  avtFVCOMParticleFileFormat.C
  avtFVCOM_MTSDFileFormat.C
  avtFVCOM_MTMDFileFormat.C
  avtCCSMFileFormat.C
  avtBOUTFileFormat.C
  avtNETCDFReaderBase.C
  avtBasicNETCDFReader.C
  avtCCSMReader.C
  avtFVCOMReader.C
  NETCDFFileObject.C
  )

ADD_VISIT_INTERFACE_READER(VisItNETCDF "1.0"
  VISIT_READERS
    NAME avtADAPTFileFormat         TYPE STSD
    NAME avtLODIFileFormat          TYPE MTSD
    NAME avtLODIParticleFileFormat  TYPE MTSD
    NAME avtFVCOM_STSDFileFormat    TYPE STSD
    NAME avtFVCOM_MTMDFileFormat    TYPE MTMD
    NAME avtFVCOM_MTSDFileFormat    TYPE MTSD
    NAME avtFVCOMParticleFileFormat TYPE MTSD
    NAME avtCCSM_MTSD_FileFormat    TYPE MTSD INCLUDE avtCCSMFileFormat
    NAME avtCCSM_STSD_FileFormat    TYPE STSD INCLUDE avtCCSMFileFormat
    NAME avtBOUTFileFormat          TYPE MTM
  VISIT_INTERFACE_CALL ${INTERFACE_MACRO_CALL}
  VISIT_INTERFACE_FILE ${INTERFACE_MACRO_FILE}
  SERVER_SOURCES ${SOURCES}
  )
