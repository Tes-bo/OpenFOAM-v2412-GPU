# petsc4Foam

## General Description

Library for [OpenFOAM] that provides a solver for embedding [Petsc][Petsc]
and its external dependencies (i.e. [Hypre][Hypre]) into arbitrary OpenFOAM simulations.

The library provides:
- ldu2csr matrix conversion.
- selection of solvers and preconditioners available in [Petsc]
  and in its external libraries (i.e. [Hypre]).
- editing of run-time options by a dictionary and/or a rc file.
- basic caching implementation of a Petsc matrix.


## License

The source code license: GPL-3.0-or-later


## Requirements

1. [OpenFOAM-v1912] or newer, or a recent [development version][OpenFOAM-git]
   from [OpenFOAM.com][OpenFOAM].
2. [Petsc] 3.10 or newer,
   optionally compiled with [Hypre] support (optionally with MPI support).


## Building

Ensure that the OpenFOAM environment is active and that Petsc
can be found (check that the `PETSC_ARCH_PATH` environment variable is properly set).
Read [How to build Petsc for OpenFOAM][Petsc-Installation] for more information
or visit the Petsc website [How to build Petsc][Petsc-Ext-Installation].


### To install in the normal OpenFOAM directories (using `wmake`)

Using the supplied `Allwmake` script without arguments:
```
./Allwmake
```
will install the library under `FOAM_MODULE_LIBBIN`.

If such variable is not defined, the standard `FOAM_USER_LIBBIN` will be used.

To install into different locations, you can use the `-prefix=PATH` or
`-prefix=shortcut` option (for OpenFOAM-v2006 and newer):
```
./Allwmake -prefix=/install/path

# Installing into the OpenFOAM FOAM_LIBBIN
./Allwmake -prefix=openfoam
```
or specify via the environment.
```
FOAM_MODULE_LIBBIN=/install/path  ./Allwmake
```

## How to use it

In order to use the library, two changes are required:
- add the libpetscFoam library to the optional keyword entry libs of the control dict file

      libs       (petscFoam);

- set the name of the solver and preconditioner in each solver of the fvSolution to petsc.
  The options database keys of each PETSc object have to be added in the petsc/options subdict of each solver equation.
  The default behaviour of the library is to convert the matrix from LDU to CSR at each time step.
  However, the user can change the cache update frequency among the following choices:

  - never (none)
  - always
  - periodic
  - adaptive

  The cache update frequency is set for both matrix and preconditioner in the petsc/caching subdict.

An example is reported below. Other examples can be found in the tutorial folder or in the [HPC repo](https://develop.openfoam.com/committees/hpc).
For more details, the user can read the paper [1].

    solvers
    {
        p
        {
            solver          petsc;
            preconditioner  petsc;

            petsc
            {
                options
                {
                    ksp_type    cg;
                    pc_type     bjacobi;
                    sub_pc_type icc;
                }

                caching
                {
                    matrix
                    {
                        update always;
                    }

                    preconditioner
                    {
                        update always;
                    }
                }
            }
        }
    }


## Authors / Contributors

| Name | Affiliation | Email
|------|-------|-----------|
| Mark Olesen  | ESI-OpenCFD | |
| Simone Bna   | CINECA | simone.bna@cineca.it |
| Stefano Zampini | KAUST | stefano.zampini@gmail.com |


----

[OpenFOAM]: https://www.openfoam.com
[OpenFOAM-v1906]: https://www.openfoam.com/releases/openfoam-v1906/
[OpenFOAM-git]: https://develop.openfoam.com/Development/openfoam

[Hypre]: https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods
[Petsc]: https://www.mcs.anl.gov/petsc/
[Petsc-Installation]: doc/README.md
[Petsc-Ext-Installation]: https://www.mcs.anl.gov/petsc/documentation/installation.html


## References

 [1] S. Bnà, I. Spisso, M. Olesen, G. Rossi *PETSc4FOAM: A Library to plug-in PETSc into the OpenFOAM
Framework* [PRACE White paper](https://prace-ri.eu/wp-content/uploads/WP294-PETSc4FOAM-A-Library-to-plug-in-PETSc-into-the-OpenFOAM-Framework.pdf)
