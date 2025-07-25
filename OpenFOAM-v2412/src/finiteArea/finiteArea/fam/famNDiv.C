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

#include "famNDiv.H"
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
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    return fa::convectionScheme<Type>::New
    (
        vf.mesh(),
        flux,
        vf.mesh().divScheme(name)
    ).ref().famDiv(flux, vf);//TODO calculate normal
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf,
    const word& name
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf, name));
    tflux.clear();

    return Div;
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const edgeScalarField& flux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fam::ndiv(flux, vf, "div("+flux.name()+','+vf.name()+')');
}


template<class Type>
tmp<faMatrix<Type>>
ndiv
(
    const tmp<edgeScalarField>& tflux,
    GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> Div(fam::ndiv(tflux(), vf));
    tflux.clear();

    return Div;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
