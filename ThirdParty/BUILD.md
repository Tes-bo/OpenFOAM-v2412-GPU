# OpenFOAM&reg; ThirdParty Build

OpenFOAM depends to a certain extent on third-party libraries
(*opensource only*). It also provides some interfaces to *opensource* or
*proprietary* libraries. This third-party collection contains configurations and
scripts for building third-party packages. It will normally only be used in
conjunction with the corresponding OpenFOAM version.

## Directory Organization

The ThirdParty directory contains a number of build scripts as well as
some directories:

| Directory         | Contains
|-------------------|-------------------------------------------------------
| etc/              | auxiliary tools and content used for the build process
| build/            | intermediate build objects
| platforms/        | the installation directories


## Configuration of Third-Party Versions

For most of the build scripts, the default software version
is provided by an appropriate OpenFOAM `etc/config.sh/...` entry.
This approach avoids duplicate entries for the default versions and
ensures the best overall consistency between the OpenFOAM installation
and its corresponding third-party installation.

Nonetheless, the distributed make scripts can generally be used for a
variety of versions of the third-party libraries, with the software
version specified on the command-line. For example,
```
$ ./makeFFTW -help
Usage: makeFFTW [OPTION] [fftw-VERSION]
```

---

## Before Starting

0. Review the [system requirements][link openfoam-require]
   and decide on the following:
   * compiler type/version - if the system compiler is not relatively recent,
     you will need a [third-party compiler](#makeGcc) installation.
   * MPI type/version.
   * ParaView type/version.
   * CMake type/version, ...
1. If you are using a system MPI (eg, openmpi), ensure that this environment
   has also been properly activated for your user.
   Often (but not always) a `mpi-selector` command is available for this purpose.
   You may need to open a new shell afterwards for the change to take effect.
   Using the following command may help diagnosing things:
```
which mpicc
```
2. Adjust the OpenFOAM `etc/bashrc`, `etc/config.sh/...` or equivalent
   the `etc/prefs.{csh,sh}` files to reflect your preferred configuration.
   For many config files, there are several configuration possibilities:
   - Define a particular third-party version.
   - Use a system installation.
   - Disable use of an optional component.
   - Define an alternative site-wide central location.

See the [OpenFOAM configuration information][link openfoam-config]
for more details.
After making the desired changes, use the `wmRefresh` alias
or equivalent to enable the configuration.


---

## Building

Many components of ThirdParty are *optional* or are invoked
automatically as part of the top-level OpenFOAM `Allwmake`.
Nonetheless it may be necessary or useful to build various
ThirdParty components prior to building OpenFOAM itself.

### Bootstrapping *(optional)*

* `makeGcc` _or_ `makeLLVM` <a name="makeGcc"></a>
   - Makes a third-party [gcc](#gcc-compiler) or [clang](#clang-compiler) installation,
     which is needed if the system gcc is [too old](#gcc-compiler).
     If your system compiler is recent enough, you can skip this step.
   - If you do use this option, you will need the following adjustments to the
     OpenFOAM `etc/bashrc` or your equivalent `prefs.sh` file:
     - `WM_COMPILER_TYPE=ThirdParty`
     - `WM_COMPILER=Gcc48` (for example)
     - `WM_COMPILER=Clang40` (for example)
     - or `WM_COMPILER=Clang` and adjust `clang_version` in the OpenFOAM
     `etc/config.sh/compiler` or equivalent.
   - More description is contained in the header comments of the
     `makeGcc` and `makeLLVM` files.
   - *Attention*: If you are building a newer version of clang, you may need to
     update your CMake beforehand.
* `makeCmake`
   - Makes a third-party [CMake](#general-packages) installation, which is
     needed if a system CMake does not exist or is [too old](#min-cmake),
   - Note that CMake is being used by an number of third-party packages
     (CGAL, LLVM, ParaView, VTK, ...)
     so this may become an increasingly important aspect of the build.

Note that the order of the bootstrapping process may need to be
reversed, or even require a few loops. For example, if you may need a
newer version of CMake before being able to build LLVM/Clang and
subsequently use the newly build clang to create a newer version of
CMake in the desired location.

Additionally, if you are using clang but with ThirdParty locations for
gmp/mpfr you will need some extra work. Here is an example:

* Compile a new ThirdParty clang version:
```
./makeLLVM llvm-4.0.1
```
* Now adjust the OpenFOAM `prefs.sh` to use the new compiler settings,
  and update the OpenFOAM environment (eg, `wmRefresh`)

* Next use (abuse) the `makeGcc` script to compile gmp/mpfr libraries.
  It is best to pass the desired versions explicitly, and necessary
  to set the CC/CXX variables so that the correct compiler is used:
```
CC=clang CXX=clang++  ./makeGcc gmp-6.1.2 mpfr-4.0.0 gcc-system
```
  specifying `gcc-system` effectively disables building of gcc,
  but will build the gmp/mpfr components.

* As a final step, it will be necessary to add the ThirdParty
  gmp/mpfr locations in the OpenFOAM config files since they are
  normally only used in combination with a ThirdParty gcc.
  The location to make these changes is in the `etc/config.sh/CGAL`,
  since this is the component that uses the mpfr library.
  For example,
```
gmp_version=gmp-6.1.2
mpfr_version=mpfr-4.0.0
export GMP_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$gmp_version
export GMP_ARCH_PATH=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER_ARCH/$mpfr_version
```
* Update update the OpenFOAM environment (eg, `wmRefresh`) again.


### Build Sequence

1. `Allwmake`
   - This will be automatically invoked by the top-level OpenFOAM `Allwmake`, but
     can also be invoked directly to find possible build errors.
   - Builds an mpi library (openmpi or mpich), scotch decomposition, boost, CGAL, FFTW.
   - If the optional kahip or metis  directories are found, they will also be compiled.
2. `makeParaView`  *(optional)*
   - Can be useful for visualization and for run-time post-processing function objects.
     You can build this at a later point in time, and rebuild the visualization modules.
3. Make any additional optional components


#### Optional Components

`makeAdios2`
- Only required for [ADIOS](#parallel) support.

`makeCGAL`
- Builds [boost](#general-packages) and [CGAL](#general-packages).
  Automatically invoked from the ThirdParty `Allwmake`,
  but can be invoked directly to resolve possible build errors.

`makeFFTW`
- Builds [FFTW](#general-packages).
  Automatically invoked from the ThirdParty `Allwmake`,
  but can be invoked directly to resolve possible build errors.

`makeKAHIP`
- Builds [KaHIP](#parallel) decomposition library.
  Automatically invoked from the ThirdParty `Allwmake`,
  but can be invoked directly to resolve possible build errors.

`makeOPENMPI`
- Builds [OPENMPI](#parallel) library.

`makePETSC`
- Only required for [PETSC](#general-packages) support.

`makeMETIS`
- Builds [METIS](#parallel) decomposition library.
  Automatically invoked from the ThirdParty `Allwmake`,
  but can be invoked directly to resolve possible build errors.

`makeSCOTCH`
- Builds [scotch, ptscotch](#parallel) decomposition libraries.
  Automatically invoked from the ThirdParty `Allwmake`,
  but can be invoked directly for special configurations
  or to resolve possible build errors.

`makeMGridGen`
- Optional agglomeration routines.

`makeCCMIO`
- Only required for conversion to/from STARCD/STARCCM+ files.

`makeMesa`, `makeVTK`
- Additional support for building offscreen rendering components.
  Useful if you want to render on computer servers without graphics cards.
  The `makeParaView.example` and `makeVTK.example` files offer some
  suggestions about compiling such a configuration.

`makeQt`
- Script to build a [Qt](#makeQt), including qmake.
- Possibly needed for `makeParaView`.
- The associated `etc/relocateQt` may be of independent use.
  Read the file for more details.

`makeGperftools`
- Build gperftools (originally Google Performance Tools)

`minCmake`
- Scour specified directories for CMakeLists.txt and their cmake_minimum.
  Report in sorted order.

`Allclean`
- After building, this script may be used to remove intermediate build
  information and save some disk space.


## Build Notes

### CGAL
- The [zlib][page zlib] library and development headers are required.

### Scotch
- The [zlib][page zlib] library and development headers are required.


### Mesa
- Needed for off-screen rendering.
- Building with [mesa-11][link mesa11] and [mesa-13][link mesa13] both
  seem okay, as does building with [mesa-17][link mesa].
- Building with mesa-12 is not possible since it fails to create
  the necessary `include/GL` directory and `osmesa.h` file.

### VTK
- Needed for off-screen rendering and run-time post-processing without
  ParaView.
- Rather than downloading VTK separately, it is easy to reuse the VTK
  sources that are bundled with ParaView.
  For example, by using a symbolic link:
```
ln -s ParaView-v5.6.0/VTK VTK-8.2.0
```
  The appropriate VTK version number can be found from the contents of
  the `vtkVersion.cmake` file.
  For example,
```
$ cat ParaView-v5.6.0/VTK/CMake/vtkVersion.cmake
```
  contains this type of information
```
# VTK version number components.
set(VTK_MAJOR_VERSION 8)
set(VTK_MINOR_VERSION 2)
set(VTK_BUILD_VERSION 0)
```

### ParaView
- Building ParaView requires CMake, qmake and a `qt` development files.
  Use the `-cmake`, `-qmake` and `-qt-*` options for `makeParaView` as
  required.
  See additional notes below about [making Qt](#makeQt) if necessary.

**NOTE** this step may not be entirely successful for your particular
system. Building ParaView itself is generally not a significant problem
but its dependency on particular Qt versions can be a problem.

If you fail at this step due to Qt dependencies, you may have success
building a slightly older Qt version. Another alternative may be to
use the ParaView-5.4 sources that were included in the ThirdParty-v1712
source pack. This, however, has never been tested in combination with
the OpenFOAM Catalyst insitu visualization.

If you encounter problems starting paraview (eg, segmentation fault)
it could be related to the graphics libraries not being properly found.
The `LIBGL_DEBUG` environment variable can be used to obtain more
detailed information:
```
$ export LIBGL_DEBUG=verbose
```

The call to paraview now emits this information:
```
libGL: screen 0 does not appear to be DRI2 capable
libGL: OpenDriver: trying /usr/lib64/dri/tls/swrast_dri.so
libGL: OpenDriver: trying /usr/lib64/dri/swrast_dri.so
...
libGL error: No matching fbConfigs or visuals found
libGL error: failed to load driver: swrast
```
Indicating that the drivers are probably not correctly installed.
The `/var/log/Xorg.0.log` file or the `dmesg` command may yield more
information to the cause.
In this particular case, the Nvidia drivers had a client/server
version mismatch.


#### 5.6.x, 5.5.x binary packages

For general functionality, the paraview version distributed with
the operating system or a [binary package][download ParaView]
may be sufficient for your needs.
- No known issues with the native OpenFOAM reader.


Using a binary package does mean you miss these elements (which may or
may not be important for you):
- cannot visualize a `blockMeshDict`
- cannot build Catalyst insitu visualization
- no alternative OpenFOAM reader module for some special features not
  found in the native OpenFOAM reader.

#### 5.5.x
- Requires patching for the vtk-m configuration (on some systems).
- Recommended patching for ParaView Catalyst.
- Recommended patching for file series
  No known issues with the native OpenFOAM reader.

#### 5.4.x
- Compiles without patching.
  No known issues with the native OpenFOAM reader.

#### 5.3.0 and older are neither recommended nor supported
- Various compilation issues and known bugs.

### Making Qt <a name="makeQt"></a>
- Building a third-party Qt installation (prior to building ParaView) requires
  some additional effort, but should nonetheless work smoothly.

1. Download a [*qt-everywhere-opensource-src*][link Qt5] package and
   unpack in the third-party directory.
2. Use the `makeQt` script with the QT version number. For example,
```
./makeQt 5.9.3
```
3. Build ParaView using this third-party QT. For example,
```
./makeParaView -qt-5.9.3  5.6.0
```
- ParaView versions prior to 5.3.0 do not properly support QT5.

- If you relocate the third-party directory to another location
  (eg, you built in your home directory, but want to install it in a
  central location), you will need to use the `etc/relocateQt` script
  afterwards.

*Note* On some older systems it can be quite difficult to build the
latest QT. In these cases, it is sometimes possible to build a
slightly older QT (eg, [qt-5.6.3][link Qt56]) instead.


---

## Versions

### Gcc Compiler <a name="gcc-compiler"></a>

The minimum version of gcc required is **4.8.5**

| Name              | Location
|-------------------|--------------------------------------------
| [gcc][page gcc]   | [releases][link gcc]
| [gmp][page gmp]   | system is often ok, otherwise [download][link gmp]
| [mpfr][page mpfr] | system is often ok, otherwise [download][link mpfr]
| [mpc][page mpc]   | system is often ok, otherwise [download][link mpc]


#### Potential MPFR conflicts

If you elect to use a third-party version of mpfr, you may experience
conflicts with your installed system mpfr.
On some systems, mpfr is compiled as *non-threaded*, whereas the
third-party will use *threaded* by default.
This can cause some confusion at the linker stage, since it may
resolve the system mpfr first (and find that it is *non-threaded*).

You can avoid this by one of two means:
1. Use system components for gmp/mpfr/mpc:  `makeGcc -system ...`
2. Use third-party mpfr, but without threading: `makeGcc -no-threadsafe ...`


#### 32-bit build (on 64-bit)

If you have a 64-bit system, but wish to have a 32-bit compiler, you
will need to enable multi-lib support for Gcc: `makeGcc -multilib`,
which is normally disabled, since many (most?) 64-bit systems do not
install the 32-bit development libraries by default.


### Clang Compiler <a name="clang-compiler"></a>

The minimum version of clang required is **3.7**

*Attention*: If you are building a newer version of clang, you may need to
update your CMake beforehand since GNU *configure* can only be used prior
to clang version 3.9.

If your system gcc is particularly old
(see [minimum gcc requirements for clang](#min-gcc))
you may have additional hurdles to using the newest versions of clang.


| Name                  | Location
|-----------------------|------------------------
| [clang][page clang]   | [download][link clang] or [newer][newer clang]
| [llvm][page llvm]     | [download][link llvm] or [newer][newer llvm]
| [openmp][page omp]    | [download][link omp] or [newer][newer omp]


### General <a name="general-packages"></a>

| Name                  | Location
|-----------------------|------------------------
| [CMake][page cmake]   | [download][link cmake]
| [boost][page boost]   | [download][link boost]
| [CGAL][page CGAL]     | [download][link CGAL]
| [FFTW][page FFTW]     | [download][link FFTW]
| [HDF5][page HDF5]     | [download][link HDF5]
| [PETSC][page PETSC]   | [download][link PETSC]
| [HYPRE][page HYPRE]   | [repo][repo HYPRE] or [download][link HYPRE]
| [ADF/CGNS][page CGNS], ccm | [link ccmio][link ccmio]
| gperftools            | [repo][repo gperftools] or [download][link gperftools]


### Parallel Processing <a name="parallel"></a>

| Name                  | Location
|-----------------------|------------------------
| [openmpi][page openmpi] | [download][link openmpi]. ***Some openmpi2/openmpi3 versions exhibit [stability issues](https://github.com/open-mpi/ompi/issues/5375)***
| [adios][page adios]   | [repo][repo adios] or [github download][link adios]
| [scotch, ptscotch][page scotch] | [repo][repo scotch] or [download][link scotch] or [older][older scotch]
| [kahip][page kahip] | [download][link kahip] or [older][older kahip]
| [metis][page metis] | [download][link metis]


### Visualization <a name="viz-version"></a>

| Name                  | Location
|-----------------------|------------------------
| [MESA][page mesa]     | [download][link mesa] or [older 13][link mesa13], [older 11][link mesa11]
| [ParaView][page ParaView] | [download][link ParaView] or older [paraview-56][link ParaView56] or [binaries][download ParaView]
| [Qt][page Qt]         | [QT5][link Qt5] for ParaView-5.3.0 and later, or the [older qt-56][link Qt56] for older systems.


### CMake Minimum Requirements <a name="min-cmake"></a>

The minimum CMake requirements for building various components.
```
2.8         llvm-3.4.2
2.8.11      CGAL-4.9
2.8.11      CGAL-4.11
3.1         CGAL-4.12 -> CGAL-4.14.3
2.8.12.2    llvm-3.7.0
2.8.12.2    llvm-3.8.0
2.8.4       cmake-3.6.0
3.3         ParaView-5.6.3
3.10        ParaView-5.7.0 - ParaView-5.8.0
3.12 - 3.16 ParaView-5.9.0 - ParaView-5.10.1
3.4.3       llvm-3.9.1
3.4.3       llvm-4.0.0 - llvm-6.0.0
3.12        ADIOS2
```

### GCC Minimum Requirements <a name="min-gcc"></a>

The minimum gcc/g++ requirements for building various components.
```
4.7         llvm-3.7.0
4.7         llvm-3.6.2
4.7         llvm-3.5.2
4.4         llvm-3.4.2
```

If your system gcc/g++ is too old to build the desired llvm/clang
version, you may need to build a lower llvm/clang version and then use
that clang compiler for building the newer llvm/clang version.

<!-- OpenFOAM, Internal links -->

[link openfoam-readme]: https://develop.openfoam.com/Development/openfoam/blob/develop/README.md
[link openfoam-config]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Config.md
[link openfoam-build]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Build.md
[link openfoam-require]: https://develop.openfoam.com/Development/openfoam/blob/develop/doc/Requirements.md
[link third-readme]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/README.md
[link third-build]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/BUILD.md
[link third-require]: https://develop.openfoam.com/Development/ThirdParty-common/blob/develop/Requirements.md


## Software Links

If you have access to this information file in its raw form, you will
note that most of the software pages and download links are listed as
a collection at the end of the file. This is designed to make it
easier to use `grep` and find the relevant pages and links.


<!-- gcc-related -->
[page gcc]:       http://gcc.gnu.org/releases.html
[page gmp]:       http://gmplib.org/
[page mpfr]:      http://www.mpfr.org/
[page mpc]:       http://www.multiprecision.org/

[link gcc]:       http://gcc.gnu.org/releases.html
[link gmp]:       ftp://ftp.gnu.org/gnu/gmp/gmp-6.2.0.tar.xz
[link mpfr]:      ftp://ftp.gnu.org/gnu/mpfr/mpfr-4.0.2.tar.xz
[link mpc]:       ftp://ftp.gnu.org/gnu/mpc/mpc-1.1.0.tar.gz


<!-- clang-related -->
[page llvm]:      http://llvm.org/
[page clang]:     http://clang.llvm.org/
[page omp]:       http://openmp.llvm.org/

[link clang]:     http://llvm.org/releases/3.7.1/cfe-3.7.1.src.tar.xz
[link llvm]:      http://llvm.org/releases/3.7.1/llvm-3.7.1.src.tar.xz
[link omp]:       http://llvm.org/releases/3.7.1/openmp-3.7.1.src.tar.xz

[newer clang]:    http://llvm.org/releases/4.0.1/cfe-4.0.1.src.tar.xz
[newer llvm]:     http://llvm.org/releases/4.0.1/llvm-4.0.1.src.tar.xz
[newer omp]:      http://llvm.org/releases/4.0.1/openmp-4.0.1.src.tar.xz


<!-- parallel -->
[page adios]:     https://csmd.ornl.gov/software/adios2
[repo adios]:     https://github.com/ornladios/ADIOS2
[link adios]:     https://github.com/ornladios/ADIOS2/archive/v2.10.1.tar.gz

[page hdf5]:      https://www.hdfgroup.org
[repo hdf5]:      https://github.com/HDFGroup/hdf5
[link hdf5]:      https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.4.3.tar.gz

[page zfp]:       http://computation.llnl.gov/projects/floating-point-compression/zfp-versions

[page scotch]:    https://www.labri.fr/perso/pelegrin/scotch/
[repo scotch]:    https://gitlab.inria.fr/scotch/scotch
[link scotch]:    https://gitlab.inria.fr/scotch/scotch/-/archive/v6.1.3/scotch-v6.1.3.tar.gz
[link scotch60_10]:  https://gforge.inria.fr/frs/download.php/file/38350/scotch_6.0.10.tar.gz
[link scotch60_9]:   https://gforge.inria.fr/frs/download.php/file/38187/scotch_6.0.9.tar.gz
[link scotch61_0]:   https://gforge.inria.fr/frs/download.php/file/38352/scotch_6.1.0.tar.gz
[link scotch61_3]:   https://gitlab.inria.fr/scotch/scotch/-/archive/v6.1.3/scotch-v6.1.3.tar.gz
[older scotch]:   https://gforge.inria.fr/frs/download.php/file/38352/scotch_6.1.0.tar.gz
[newer scotch]:   https://gitlab.inria.fr/scotch/scotch/-/archive/v7.0.1/scotch-v7.0.1.tar.gz

[page kahip]:     http://algo2.iti.kit.edu/documents/kahip/
[older kahip]:    https://github.com/KaHIP/KaHIP/archive/v3.14.tar.gz
[link kahip]:     https://github.com/KaHIP/KaHIP/archive/v3.15.tar.gz

[page metis]:     http://glaros.dtc.umn.edu/gkhome/metis/metis/overview
[link metis]:     http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz

[page openmpi]:   http://www.open-mpi.org/
[older openmpi]:  https://download.open-mpi.org/release/open-mpi/v1.10/openmpi-1.10.7.tar.bz2
[link openmpi]:   https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.2.tar.bz2

[page mpich]:     http://www.mpich.org/
[link mpich]:     http://www.mpich.org/static/downloads/3.4.2/mpich-3.4.2.tar.gz

[page mvpapich]:  http://mvapich.cse.ohio-state.edu/
[link mvpapich]:  https://mvapich.cse.ohio-state.edu/download/mvapich/mv2/mvapich2-2.3.7-1.tar.gz

<!-- general -->
[page cmake]:     http://www.cmake.org/
[link cmake]:     https://cmake.org/files/v3.8/cmake-3.8.2.tar.gz

[page boost]:     http://boost.org
[link boost]:     https://sourceforge.net/projects/boost/files/boost/1.74.0/boost_1_74_0.tar.bz2

[page CGAL]:      http://cgal.org
[link CGAL]:      https://github.com/CGAL/cgal/releases/download/releases/CGAL-4.14.3/CGAL-4.14.3.tar.xz

[page FFTW]:      http://www.fftw.org/
[link FFTW]:      http://www.fftw.org/fftw-3.3.10.tar.gz

[page petsc]:     https://petsc.org/
[link petsc]:     https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-lite-3.21.2.tar.gz

[page hypre]:     https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods/
[repo hypre]:     https://github.com/hypre-space/hypre/
[link hypre]:     https://github.com/hypre-space/hypre/archive/v2.31.0.tar.gz

[page cgns]:      http://cgns.github.io/
[link ccmio]:     http://portal.nersc.gov/project/visit/third_party/libccmio-2.6.1.tar.gz (check usage conditions)
[altlink ccmio]:  https://sourceforge.net/projects/foam-extend/files/ThirdParty/libccmio-2.6.1.tar.gz (check usage conditions)

[repo gperftools]: https://github.com/gperftools/gperftools
[link gperftools]: https://github.com/gperftools/gperftools/releases/download/gperftools-2.5/gperftools-2.5.tar.gz

[page zlib]:      https://www.zlib.net/
[link zlib]:      https://sourceforge.net/projects/libpng/files/zlib/1.2.11/zlib-1.2.11.tar.xz


<!-- profilint -->

[link extrae]:     https://github.com/bsc-performance-tools/extrae/archive/refs/tags/4.1.7.tar.gz
[link libunwind]:  https://github.com/libunwind/libunwind/releases/download/v1.8.1/libunwind-1.8.1.tar.gz
[link papi]:       http://icl.utk.edu/projects/papi/downloads/papi-7.0.1.tar.gz

<!-- Visualization -->

[page ParaView]:  http://www.paraview.org/
[download ParaView]: https://www.paraview.org/download/
[link ParaView56]: http://www.paraview.org/files/v5.6/ParaView-v5.6.3.tar.xz
[link ParaView57]: http://www.paraview.org/files/v5.7/ParaView-v5.7.0.tar.xz
[link ParaView58]: http://www.paraview.org/files/v5.8/ParaView-v5.8.1.tar.xz
[link ParaView59]: http://www.paraview.org/files/v5.9/ParaView-v5.9.1.tar.xz
[link ParaView510]: http://www.paraview.org/files/v5.10/ParaView-v5.10.1.tar.xz
[link ParaView511]: http://www.paraview.org/files/v5.11/ParaView-v5.11.2.tar.xz
[link ParaView512]: http://www.paraview.org/files/v5.12/ParaView-v5.12.1.tar.xz
[link ParaView]:    http://www.paraview.org/files/v5.12/ParaView-v5.12.1.tar.xz

[page mesa]:  http://mesa3d.org/
[link mesa]:  ftp://ftp.freedesktop.org/pub/mesa/mesa-17.1.1.tar.xz
[link mesa13]: ftp://ftp.freedesktop.org/pub/mesa/13.0.6/mesa-13.0.6.tar.xz
[link mesa11]: ftp://ftp.freedesktop.org/pub/mesa/older-versions/11.x/11.2.2/mesa-11.2.2.tar.xz

[page Qt]: https://www.qt.io/download-open-source/
[repo Qt]: http://code.qt.io/cgit/qt-creator/qt-creator.git
[link Qt56]: http://download.qt.io/official_releases/qt/5.6/5.6.3/single/qt-everywhere-opensource-src-5.6.3.tar.xz
[link Qt5]: http://download.qt.io/official_releases/qt/5.9/5.9.3/single/qt-everywhere-opensource-src-5.9.3.tar.xz

---

Copyright 2016-2024 OpenCFD Ltd
