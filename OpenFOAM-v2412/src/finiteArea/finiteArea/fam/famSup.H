/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Foam::fam

Description
    Calculate the finiteArea matrix for implicit and explicit sources.

SourceFiles
    famSup.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_famSup_H
#define Foam_famSup_H

#include "areaFieldsFwd.H"
#include "faMatrix.H"
#include "zeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Namespace fam Functions
\*---------------------------------------------------------------------------*/

namespace fam
{
   // Explicit source

        //- A no-op source
        template<class Type>
        zeroField Su
        (
            const Foam::zero,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<faMatrix<Type>> Su
        (
            const dimensioned<Type>& su,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Su
        (
            const DimensionedField<Type, areaMesh>& su,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Su
        (
            const tmp<DimensionedField<Type, areaMesh>>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Su
        (
            const tmp<GeometricField<Type, faPatchField, areaMesh>>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );


    // Implicit source

        //- A no-op source
        template<class Type>
        zeroField Sp
        (
            const Foam::zero,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<faMatrix<Type>> Sp
        (
            const dimensionedScalar& sp,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Sp
        (
            const DimensionedField<scalar, areaMesh>& sp,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Sp
        (
            const tmp<DimensionedField<scalar, areaMesh>>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> Sp
        (
            const tmp<GeometricField<scalar, faPatchField, areaMesh>>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );


    // Implicit/Explicit source depending on sign of coefficient

        //- A no-op source
        template<class Type>
        zeroField SuSp
        (
            const Foam::zero,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        //- A uniform source (no-op for small values)
        template<class Type>
        tmp<faMatrix<Type>> SuSp
        (
            const dimensionedScalar& susp,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> SuSp
        (
            const DimensionedField<scalar, areaMesh>& susp,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> SuSp
        (
            const tmp<DimensionedField<scalar, areaMesh>>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );

        template<class Type>
        tmp<faMatrix<Type>> SuSp
        (
            const tmp<areaScalarField>&,
            const GeometricField<Type, faPatchField, areaMesh>&
        );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "famSup.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
