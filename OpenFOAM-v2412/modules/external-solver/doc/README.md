# PETSC

Configure and install PETSc.
For example,
```
./configure
    --with-64-bit-indices=0 \
    --with-precision=double \
    --prefix=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_LABEL_OPTION/petsc-git \
    PETSC_ARCH=$WM_OPTIONS

make PETSC_DIR=... PETSC_ARCH=... all
make PETSC_DIR=... PETSC_ARCH=... install
```

Note that this approach does not account for mixing of mpi
distributions or versions. Also, it will configure PETSc in debug mode, with an impact on performances.

Export the appropriate PETSC_ARCH_PATH location and/or adjust config
files(s):

- global: `$WM_PROJECT_DIR/etc/config.sh/petsc`
- user:   `$HOME/.OpenFOAM/config.sh/petsc`



Another possible configuration for an optimized build of PETSc
```
./configure
    --with-64-bit-indices=0 \
    --with-precision=double \
    --with-debugging=0 \
    --COPTFLAGS=-O3 \
    --CXXOPTFLAGS=-O3 \
    --FOPTFLAGS=-O3 \
    --prefix=$WM_THIRD_PARTY_DIR/platforms/$WM_ARCH$WM_COMPILER$WM_PRECISION_OPTION$WM_LABEL_OPTION/petsc-git \
    PETSC_ARCH=$WM_OPTIONS \
    --download-openblas \
    --with-fc=0 \
    --with-mpi-dir=$MPI_ARCH_PATH
```

The ThirdParty [makePETSC][makePETSC] script may be of some use.

For performant runs, we strongly advise you to specify the location of your favorite BLAS/LAPACK libraries either using
```
--with-blaslapack-lib=_comma_separated_list_of_libraries_needed_to_resolve_blas/lapack_symbols
```
or
```
--with-blaslapack-dir=_prefix_location_where_to_find_include_and_libraries
```

For additional information on the PETSc configure process, see [here](https://www.mcs.anl.gov/petsc/documentation/installation.html).

# Running

Requires LD_LIBRARY_PATH to include PETSc information.
For example, using the config settings:
```
eval $(foamEtcFile -sh -config petsc -- -force)
```

And an additional library loading in the `system/controlDict`.
For example,
```
libs    (petscFoam);
```
Alternatively, can preload this library directly from the command-lib
for the solver.
For example,
```
simpleFoam -lib petscFoam ...
```

----

[makePETSC]: https://develop.openfoam.com/Development/ThirdParty-common/-/blob/master/makePETSC
