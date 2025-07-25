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

#include "steadyStateFaDdtScheme.H"
#include "facDiv.H"
#include "faMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensioned<Type> dt
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt("+dt.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        dt.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensioned<Type> dt
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt("+dt.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        dt.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt("+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        vf.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt0("+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        vf.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt("+rho.name()+','+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        rho.dimensions()*vf.dimensions()/dimTime
    );
}

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        rho.dimensions()*vf.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt("+rho.name()+','+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        rho.dimensions()*vf.dimensions()/dimTime
    );
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
steadyStateFaDdtScheme<Type>::facDdt0
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return GeometricField<Type, faPatchField, areaMesh>::New
    (
        "ddt0("+rho.name()+','+vf.name()+')',
        IOobject::REGISTER,
        mesh(),
        Foam::zero{},
        rho.dimensions()*vf.dimensions()/dimTime
    );
}

template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}


template<class Type>
tmp<faMatrix<Type>>
steadyStateFaDdtScheme<Type>::famDdt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    tmp<faMatrix<Type>> tfam
    (
        new faMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimArea/dimTime
        )
    );

    return tfam;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
