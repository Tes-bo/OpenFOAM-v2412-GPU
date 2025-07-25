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

\*---------------------------------------------------------------------------*/

#include "famDiv.H"
#include "faMesh.H"
#include "faMatrix.H"
#include "faConvectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<faMatrix<Type>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    const areaVectorField& n = vf.mesh().faceAreaNormals();

    tmp<faMatrix<Type>> tM
    (
        fa::convectionScheme<Type>::New
        (
            vf.mesh(),
            flux,
            vf.mesh().divScheme(name)
        ).ref().famDiv(flux, vf)
    );
    faMatrix<Type>& M = tM.ref();

    GeometricField<Type, faPatchField, areaMesh> v
    (
        fa::convectionScheme<Type>::New
        (
            vf.mesh(),
            flux,
            vf.mesh().divScheme(name)
        ).ref().facDiv(flux, vf)
    );

    //HJ Check if the product is from left or right.  HJ, 6/Dec/2016
    M -= n*(n & v);

    return tM;
}


template<class Type>
tmp<faMatrix<Type>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> Div(fam::div(tflux(), vf, name));
    tflux.clear();

    return Div;
}


template<class Type>
tmp<faMatrix<Type>>
div
(
    const edgeScalarField& flux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::div(flux, vf, "div("+flux.name()+','+vf.name()+')');
}

template<class Type>
tmp<faMatrix<Type>>
div
(
    const tmp<edgeScalarField>& tflux,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> Div(fam::div(tflux(), vf));
    tflux.clear();

    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
