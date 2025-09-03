# OpenFOAM&reg; ThirdParty

OpenFOAM depends to a certain extent on third-party libraries
(*opensource only*). It also provides some interfaces to *opensource* or
*proprietary* libraries. This third-party package contains configurations and
scripts for building third-party packages. It should normally only be used in
conjunction with the corresponding OpenFOAM version.

## Configuration of Third-Party Versions

For most of the build scripts, the default software version
is provided by an appropriate OpenFOAM `etc/config.sh/...` entry.
This approach avoids duplicate entries for the default versions and
ensures the best overall consistency between the OpenFOAM installation
and its corresponding third-party installation.

Nonethess, the distributed make scripts can generally be used for a
variety of versions of the third-party libraries, with the software
version specified on the command-line.

---

## Before Starting

0. Review the [system requirements][link openfoam-require]
   and decide on the following:
   * compiler type/version (you may need a third-party compiler installation).
   * MPI type/version.
   * ParaView type/version.
   * CMake type/version, ...
1. Adjust the OpenFOAM `etc/bashrc`, `etc/config.sh/...` or equivalent
   `prefs.sh` files to reflect your preferred [configuration][link openfoam-config].
2. Source the updated OpenFOAM environment

---

## Building

Many components of ThirdParty are *optional* or are invoked
automatically as part of the top-level OpenFOAM `Allwmake`.
Nonetheless it may be necessary or useful to build particular
ThirdParty components prior to building OpenFOAM itself.

### Sequence

1. `makeGcc` _or_ `makeLLVM` *(optional)*
2. `makeCmake`  *(optional)*
3. `Allwmake`
   - This will be automatically invoked by the top-level OpenFOAM `Allwmake`.
4. `makeParaView`  *(optional)*
5. Any other additional optional components

If the `Allwmake` is being invoked directly (not from the OpenFOAM `Allwmake`),
it may be necessary to bootstrap the binary parts of the `wmake` toolchain
manually before proceeding. This is done with the following command:
```
$WM_PROJECT_DIR/wmake/src/Allmake
```

**It normally advisable to skip the ParaView compilation initially**.

Building ParaView from source tends to be the most difficult part of
any third-party compilation.

For general functionality, the paraview version distributed with
the operating system or a [binary package][download ParaView]
may be sufficient for your needs.


### Details

- More details can be found the [ThirdParty BUILD.md][link third-build] information.
- Some configuration details can be found in the
  [OpenFOAM doc/Config.md][link openfoam-config] information.


### Other

When building other libraries directly, with autoconfig or cmake, it
will be useful or necessary to use the same compiler and compiler
settings as are used by OpenFOAM itself. These are obtained from the
`wmake` show options and can be used to set corresponding environment
variables. For example,
```
CC="$(wmake -show-c)" CFLAGS="$(wmake -show-cflags)" ./configure
```

Here is the correspondence to commonly used environment variables

| Env variable      | Obtaining from wmake      | Meaning               |
|-------------------|---------------------------|-----------------------|
| CC                | `wmake -show-c`           | C compiler            |
| CFLAGS            | `wmake -show-cflags`      | C compiler flags      |
| CXX               | `wmake -show-cxx`         | C++ compiler          |
| CXXFLAGS          | `wmake -show-cxxflags`    | C++ compiler flags    |
| | `wmake -show-cflags-arch`   | Architecture information when linking |
| | `wmake -show-cxxflags-arch` | Architecture information when linking |


In some situations it can also be useful to have the compiler and flags
together (similar to `mpicc -show` and `mpicxx -show`):
```
wmake -show-compile-c
wmake -show-compile-cxx
```

<!-- Quick links -->

[download ParaView]: https://www.paraview.org/download/

<!-- OpenFOAM -->

[link openfoam-readme]: https://develop.openfoam.com/Development/openfoam/blob/develop/README.md
[link openfoam-config]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Config.md
[link openfoam-build]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md
[link openfoam-require]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/Requirements.md

---
Copyright 2016-2024 OpenCFD Ltd
