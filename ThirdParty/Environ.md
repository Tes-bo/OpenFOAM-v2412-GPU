## OpenFOAM&reg; ThirdParty Environment Variables

When using these scripts for building third-party libraries
various environment variables are required or expected.

### Mandatory

- **WM_THIRD_PARTY_DIR** : Location of third-party library sources etc
- **WM_PROJECT_DIR** : The OpenFOAM project directory. <br>
  Used for `wmake` information, config files, etc
- **WM_ARCH** : The target architecture (eg, Linux64)
- **WM_OSTYPE** : The operating system type (eg, POSIX, MSwindows) <br>
  Used as fallback for determining the ending for dynamic libraries.
  Used for Windows-specific handling in scotch builds.
- **WM_COMPILER** : The compiler name (eg, Gcc, Clang)
- **WM_COMPILER_LIB_ARCH** : The target library ending (eg, 64). <br>
  Typically used as `lib$WM_COMPILER_LIB_ARCH` to generate the name
  `lib64`, for example.
- **WM_LABEL_SIZE** : The OpenFOAM size for `label` (eg, 32 for int32)
- **WM_PRECISION_OPTION** : The OpenFOAM representation for `scalar`
  (eg, DP for double-precision)
- **FOAM_EXT_LIBBIN** : Common _pool_ for the library targets.


### MPI-related

- **FOAM_MPI** : The OpenFOAM MPI sub-directory name (eg, sys-openmpi)
- **WM_MPLIB** : The canonical OpenFOAM MPI name (eg, SYSTEMOPENMPI)
- **MPI_ARCH_PATH** : Location of the current MPI implementation.
  Used for additional linkage, but primarily to set MPI_ROOT
  for a cmake hint.


### Other

- **WM_CONTINUE_ON_ERROR** : Ignore build errors.
  Mostly set internally from `wmake -k`.
- **WM_NCOMPPROCS** : The number of build processes to launch.
  Often set internally from `wmake -j`.
- **WM_COMPILER_TYPE** : Used internally for a few build scripts.
  Has values of `system` or `ThirdParty`.


## Output

Similar to OpenFOAM itself, the build and final output are separated
into different directories that contain a non-colliding hierarchy of
target names. However, in contrast to the OpenFOAM target names,
there are several different types of output targets:

- *ARCH* : (eg, `linux64`) <br>
  Base infrastructure such as compiler and base libraries.
- *ARCH+COMPILER* : (eg, `linux64Gcc`) <br>
  Compiler-specific libraries and applications.
- *ARCH+COMPILER+SIZING* : (eg, `linux64GccDPInt32`) <br>
  Libraries and applications specific to a particular combination
  of OpenFOAM scalar/label sizes. This is the same as the
  ***WM_OPTIONS*** value.

By default, the locations of the build and final targets are hard-coded
relative to the WM_THIRD_PARTY_DIR:

- WM_THIRD_PARTY_DIR<b>/build/</b> : intermediate build artifacts
- WM_THIRD_PARTY_DIR<b>/platforms/</b> : installation root

The expert user or sys-admin can influence these locations with the
following environment variables

- **FOAM_THIRD_PARTY_BUILDROOT** :
  Replaces WM_THIRD_PARTY_DIR as the root for build/ and platforms/.<br>
  ***!! The specified build-root directory must exist !!***

- **FOAM_THIRD_PARTY_SOURCES** :
  Provide alternative location to WM_THIRD_PARTY_DIR/sources
  for finding source bundles.<br>

---
