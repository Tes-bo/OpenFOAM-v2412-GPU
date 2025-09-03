## General Description

Library for [OpenFOAM] that provides a runtime-selectable
function object for embedding [ParaView Catalyst][Catalyst]
in-situ visualization into arbitrary OpenFOAM simulations.

Supports in-situ conversion of the following types:
- finite volume meshes and fields. Single or multi-region.
- finite area meshes and fields. Single region.
- lagrangian (clouds). Single or multiple clouds.

### Overset

For simulations with overset meshes, internal blanking is used to hide
the *holes* so that user visualization pipelines do not require any
thresholding workarounds. Note, however, that visualization artifices
may be present in multi-overlapping regions.


## Requirements

1. [OpenFOAM-v1806] or a recent [development version][OpenFOAM-git]
   from [OpenFOAM.com][OpenFOAM].
2. [ParaView] or ParaView [Catalyst] 5.5 or newer,
   compiled with python support (optionally with MPI support).


### Patching

It is **highly recommended** to patch the ParaView 5.5 sources to
include changes ([MR2433], [MR2436]) that will be part of the
ParaView 5.6 release.

These patches are necessary to ensure that the in-situ results are
placed in the correct output directory.
Without these patches, the results will always land in the current
working directory: not in the case-local `insitu` subdirectory
and ignoring the OpenFOAM `-case` parameter.

These [patches] will be automatically applied when the `makeParaView`
script from OpenFOAM ThirdParty is used.


## Building

Ensure that the OpenFOAM environment is active and that ParaView or Catalyst
can be found (Eg, the `ParaView_DIR` environment is properly set).


### To install in the normal OpenFOAM directories (using `wmake`)

Simply use the supplied `Allwmake` script:
````
./Allwmake
````
This will install the library under
`$WM_PROJECT_DIR/platforms/$WM_OPTIONS/lib`,
which corresponds to the `$FOAM_LIBBIN` location.


### To install in arbitrary locations (using `cmake`)

Without parameters, it installs to `/usr/local` as the default location.
This can be changed as required with the cmake `CMAKE_INSTALL_PREFIX`
parameter.

````
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=/install/path ../src/catalyst
````

The installation using the `./Allwmake` script can be replicated with
either of these commands:

````
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=$WM_PROJECT_DIR/platforms/$WM_OPTIONS ../src/catalyst
make install
````
or
````
mkdir build; cd build
cmake -DCMAKE_INSTALL_PREFIX=${FOAM_LIBBIN%/*} ../src/catalyst
make install
````


### Where to start

Look at the tutorials cases to get an idea of how different things
might be done. Know what types of inputs your particular Catalyst
function object is producing. The `writeAll.py` script included as
part of the normal OpenFOAM distribution provides a good starting point
for examing your output, but also for generating initial data sets for
your pipelines. Use the script generation capabilities within ParaView
to generate your pipelines. Save a backup of your pipelines as `pvsm`
state files so that you can go back and tweak your pipelines without
manually editing the python scripts and introducing possible errors.

Be aware that there may be some combinations of pipelines or
annotations that misbehave. If you encounter such a problem, try
backing out various pieces to isolate the cause. If it does appear to be
a reproducible error, test it on a tutorial example too.


## Authors

- Mark Olesen | <mark.olesen@esi-group.com>  | (ESI-OpenCFD)
- Simone Bna  | <simone.bna@cineca.it>       | (CINECA)


## License

Licensed under GNU General Public License <http://www.gnu.org/licenses/>
with the same terms as OpenFOAM itself.

----

[OpenFOAM]: https://www.openfoam.com
[OpenFOAM-v1806]: https://www.openfoam.com/releases/openfoam-v1806/
[OpenFOAM-git]: https://develop.openfoam.com/Development/OpenFOAM-plus
[patches]:  https://develop.openfoam.com/Development/ThirdParty-plus/raw/develop/etc/patches/paraview-5.5.2


[ParaView]: https://www.paraview.org/
[Catalyst]: https://www.paraview.org/in-situ/
[MR2433]: https://gitlab.kitware.com/paraview/paraview/merge_requests/2433
[MR2436]: https://gitlab.kitware.com/paraview/paraview/merge_requests/2436
