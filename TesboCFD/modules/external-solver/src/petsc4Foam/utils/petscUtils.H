/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2019 Simone Bna
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Namespace
    Foam::PetscUtils

Description
    Utilities for OpenFOAM / PETSc interaction

SourceFiles
    petscUtils.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscUtils_H
#define Foam_module_petscUtils_H

// OpenFOAM
#include "scalarField.H"

// PETSc
#include "petscvec.h"
#include "petscksp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//- Global average value in PETSc vector
solveScalar gAverage(Vec input);

//- Global sum of PETSc vector
solveScalar gSum(Vec input);


namespace PetscUtils
{

/*---------------------------------------------------------------------------*\
                             Namespace PetscUtils
\*---------------------------------------------------------------------------*/

//- Calculate the OpenFOAM normFactor (akin to an L1 norm)
scalar normFactor
(
    Vec AdotPsi,
    Vec psi,
    Vec source,
    Vec ArowsSum
);

PetscErrorCode foamKSPConverge
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    KSPConvergedReason* reason,
    void* cctx
);

PetscErrorCode foamKSPMonitorFoam
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    void* cctx
);

PetscErrorCode foamKSPMonitorRecordInit
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    void* cctx
);

//- Set a single PETSc flag.
//  The key already contains a leading '-'
void setFlag
(
    const word& key,
    const word& val,
    const bool verbose = false
);

//- Set PETSc flags.
//  The prefix does not contain a leading '-'
void setFlags
(
    const word& prefix,
    const dictionary& dict,
    const bool verbose = false
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PetscUtils
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
