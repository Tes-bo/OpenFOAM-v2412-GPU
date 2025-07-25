/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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

InNamespace
    Foam::fvm

Description
    Calculate the finiteVolume matrix for implicit and explicit sources.

SourceFiles
    fvmSup.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_fvmSup_H
#define Foam_fvmSup_H

#include "volFieldsFwd.H"
#include "fvMatrix.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Namespace fvm Functions
\*---------------------------------------------------------------------------*/

namespace fvm
{
   // Explicit source

        //- A no-op source
        template<class Type>
        zeroField Su
        (
            const Foam::zero,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const dimensioned<Type>& su,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const DimensionedField<Type, volMesh>& su,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<DimensionedField<Type, volMesh>>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Su
        (
            const tmp<GeometricField<Type, fvPatchField, volMesh>>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Implicit source

        //- A no-op source
        template<class Type>
        zeroField Sp
        (
            const Foam::zero,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const dimensionedScalar& sp,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const DimensionedField<scalar, volMesh>& sp,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<DimensionedField<scalar, volMesh>>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> Sp
        (
            const tmp<volScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        //- A no-op source
        template<class Type>
        zeroField SuSp
        (
            const Foam::zero,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const dimensionedScalar& susp,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const DimensionedField<scalar, volMesh>& susp,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<DimensionedField<scalar, volMesh>>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );

        template<class Type>
        tmp<fvMatrix<Type>> SuSp
        (
            const tmp<volScalarField>&,
            const GeometricField<Type, fvPatchField, volMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvmSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
