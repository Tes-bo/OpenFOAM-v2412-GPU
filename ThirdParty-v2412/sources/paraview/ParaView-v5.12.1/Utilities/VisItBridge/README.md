Introduction
===============

This article describes the Visit Database Bridge ParaView plugin. The motivation for the bridge is to allow ParaView to make use of VisIt's IO components, and to explore the re-usable capabilities of VisIt and its underlying pipeline library, avt. Like ParaView, VisIt is a scientific data visualization application based on VTK and Qt. There are a number of subtle differences in the internal workings of the two applications. However, in general there are more similarities than differences. Hence porting functionality from one to the other is fairly straight forward and potentially very beneficial as both applications have large, publicly available, thoroughly tested code bases.

VisIt IO components are called database plugins. The VisIt developers have provided a consistent interface both for database plugin via the avtFileFormat interface. The VisItBridge now interfaces with the avtFileFormat to create a VTK pipeline object that emulates the VisIt database loading process. The reason for this lower level interfacing was to reduce the number of VisIt dependencies, eliminate the need for runtime loading of readers, and have faster access to the underlying VTK data. The wrapping of each reader is now done automatically be auto generating the VTK class at configure time. This allows developers to easily add more VisItReaders or even to add VisItReaders into plugins.

The new VisItBridge secondary goal was to provide an integrated build with ParaView. This was goal was completed by forking a subset of the VisIt source code and patching it for compilation against VTK head. The switching of VisIt to CMake made this possible as it allows us to configure VisIt to target the HDF5 and XDMF being built by ParaView.

Classes
---------------

The structure of our bridge is as follows:

* vtkAvtFileFormatAlgorithm – Parent VTK algorithm that implements the VTK pipeline required of a ParaView reader. Contains an avtFileFormat, avtDatabaseMetaData, and avtVariableCache.
The following classes are subclasses that handle specific file formats:

- vtkAvtSTSDFileFormatAlgorithm
- vtkAvtSTMDFileFormatAlgorithm
- vtkAvtMTSDFileFormatAlgorithm
- vtkAvtMTMDFileFormatAlgorithm

Which of VisIt's databases are available?
---------------

This depends largely on how VisIt is built because databases only are built if their dependencies are found. The following list reflects what is available in a full VisItBridge (based on VisIT 2.1.1) build. The number of new formats is 60, with the potential to grow as we add support for more third party libraries.

* ANALYZE
* ANSYS
* AUXFile
* BOV
* CEAucd
* Chombo
* Claw
* CMAT
* CTRL
* Curve2D
* DDCMD
* Dyna3D
* Enzo
* ExtrudedVol
* FLASH
* Fluent
* Gadget
* GGCM
* GTC
* GULP
* H5Nimrod
* Image
* Lines
* M3D
* M3DC1
* Mili
* Miranda
* MM5
* NASTRAN
* Nek5000
* NETCDF
* OpenFOAM
* OVERFLOW
* ParaDIS
* PATRAN
* PFLOTRAN
* Pixie
* PlainText
* PLOT2D
* PLOT3D
* Point3D
* ProteinDataBank
* RAW
* SAMRAI
* SAR
* SAS
* Silo
* Spheral
* Tecplot
* Tetrad
* TFT
* TSurf
* UNIC
* VASP
* Velodyne
* Vs
* XDMF
* XYZ

Building ParaView
===============

To build ParaView with the VisItBridge enabled you will need to set the following CMake variables:
```
 PARAVIEW_USE_VISITBRIDGE
 Boost_INCLUDE_DIR
```
Silo Support
---------------

To enable Silo readers you will need to set:
```
 VISIT_BUILD_READER_Silo
 SILO_INCLUDE_DIR
 SILO_LIBRARY
```

Mili Support
---------------

To enable Mili readers you will need to set:
```
 VISIT_BUILD_READER_Mili
 Mili_INCLUDE_DIR
 Mili_LIBRARY
```

Custom Plugin of VisIt Reader
===============

One of the major benefits of the new VisItBridge is the ability to allow developers to create plugins of any VisIt reader for ParaView. Adding the reader is incredibly easy as all you need is a CMakeLists.txt file and the reader source.

The CMakeLists.txt looks as follows where avtMyReaderFileFormat.cxx is the source of the VisIt reader, VISIT_READER_NAME is the name of the class inside the source file, VISIT_READER_TYPE is the type of file format the reader is, and the VISIT_READER_FILE_PATTERN is a space delimitated string of file extensions that the reader supports.

```cmake
PROJECT(VisItReaderExample)

cmake_minimum_required(VERSION 2.8)

FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})

SET(SOURCES
 avtMyReaderFileFormat.C
 )

ADD_VISIT_PLUGIN_READER(VisItReaderExample "1.0"
  VISIT_READER_NAME "avtMyReaderFileFormat"
  VISIT_READER_TYPE "STMD"
  VISIT_READER_FILE_PATTERN "cas"
  SERVER_SOURCES ${SOURCES}
  )

TARGET_LINK_LIBRARIES(VisItReaderExample vtkVisItAVTAlgorithms)
```

What's Not Done
===============

* We need to map VisIt ghost cells to VTK ghost cells, not doing so is going to introduce block  artifacts.
* Implement a SIL interface in the GUI. A proper SIL needs to be added to clarify the current Material panel, and to properly show what properties each block has.
* Finish Expressions. See avtExpresssionEvaluationFilter for adding support.
* AMR data. AMR data has been partially implemented. The mapping of ghost arrays has been implemented but hasn’t been fully tested.

References
==========

[Source](http://www.paraview.org/Wiki/VisIt_Database_Bridge)

License
=======

Most of VisitBridge source code is coming from Visit,
so VisitBridge is distributed under the same licence.
See LICENSE file.
