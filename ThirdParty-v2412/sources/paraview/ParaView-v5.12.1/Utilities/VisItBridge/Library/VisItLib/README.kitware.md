# VisIt fork for ParaView

This branch contains changes required to embed VisIt into ParaView via
VisItBridge. This includes changes made primarily to the build system to allow
it to be embedded into another source tree as well as a header to facilitate
mangling of the symbols to avoid conflicts with other copies of the library
within a single process.

  * Add attributes to pass commit checks within VisItBridge.
  * Add missing newlines to the ends of files.
  * Replace non-utf8 characters in imported files.
  * Rework CMake code to build the subset required for ParaView.
  * Add CMake code to integrate with ParaView.
  * Modify readers to use VTK's third party modules.
  * Ignore `DBIO_ONLY` in the `avtMath` library.
  * Prevent duplicate symbols within the OpenFOAM reader.
  * Add methods to `avtFileFormat` for use with the ParaView glue code.
  * Have the `paraDIS_tecplot` reader share `RC_cpp_lib` with the `paraDIS`
    reader.
  * Namespace various `Atom` structures to avoid symbol collisions within
    ParaView.
  * Add the `ValidNETCDFFile.h` header.
  * Add CMake lists for databases imported into ParaView.
  * Comment out unused functions in `RC_cpp_lib`.
  * Updating avtANSYSFileFormat reader to read in 10-noded tetras.
  * Read PFLOTRAN files with big-endian floating-point arrays.
  * Port MFIXCDF from netcdfcpp to netcdf
