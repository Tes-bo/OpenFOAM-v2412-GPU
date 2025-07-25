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

#include "facDdt.H"
#include "faMesh.H"
#include "faDdtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fac
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
ddt
(
    const dimensioned<Type> dt,
    const faMesh& mesh
)
{
    return fa::faDdtScheme<Type>::New
    (
        mesh,
        mesh.ddtScheme("ddt(" + dt.name() + ')')
    ).ref().facDdt(dt);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
ddt
(
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().facDdt(vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
ddt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme
        (
            "ddt(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().facDdt(rho, vf);
}


template<class Type>
tmp<GeometricField<Type, faPatchField, areaMesh>>
ddt
(
    const areaScalarField& rho,
    const GeometricField<Type, faPatchField, areaMesh>& vf
)
{
    return fa::faDdtScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme
        (
            "ddt(" + rho.name() + ',' + vf.name() + ')'
        )
    ).ref().facDdt(rho, vf);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fac

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
