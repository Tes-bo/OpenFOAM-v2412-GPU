/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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
    Calculate the matrix for the laplacian of the field.

SourceFiles
    famLaplacian.C

\*---------------------------------------------------------------------------*/

#ifndef famLaplacian_H
#define famLaplacian_H

#include "areaFieldsFwd.H"
#include "edgeFieldsFwd.H"
#include "faMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Namespace fam functions Declaration
\*---------------------------------------------------------------------------*/

namespace fam
{
    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const dimensionedScalar&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const dimensionedScalar&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const areaScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<areaScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );
    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<areaScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const edgeScalarField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<edgeScalarField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const areaTensorField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<areaTensorField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const edgeTensorField&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<edgeTensorField>&,
        const GeometricField<Type, faPatchField, areaMesh>&
    );


    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const edgeTensorField&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );

    template<class Type>
    tmp<faMatrix<Type>> laplacian
    (
        const tmp<edgeTensorField>&,
        const GeometricField<Type, faPatchField, areaMesh>&,
        const word& name
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "famLaplacian.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
