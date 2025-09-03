/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "strainRateFunction.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(strainRateFunction, 0);

    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        strainRateFunction,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::strainRateFunction::
strainRateFunction
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    strainRateFunctionCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    strainRateFunction_
    (
        Function1<scalar>::New("function", strainRateFunctionCoeffs_)
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::
strainRateFunction::read
(
    const dictionary& viscosityProperties
)
{
    generalizedNewtonianViscosityModel::read(viscosityProperties);

    strainRateFunctionCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    strainRateFunction_.reset();
    strainRateFunction_ = Function1<scalar>::New
    (
        "function",
        strainRateFunctionCoeffs_
    );

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::strainRateFunction::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    auto tnu = volScalarField::New
    (
        IOobject::scopedName(type(), IOobject::groupName("nu", nu0.group())),
        nu0.mesh(),
        dimensionedScalar(dimViscosity, Zero)
    );

    // Set the database - could not be done in constructor
    const_cast<Function1<scalar>&>(strainRateFunction_()).resetDb
    (
        nu0.mesh().thisDb()
    );

    tnu.ref().primitiveFieldRef() = strainRateFunction_->value(strainRate);

    volScalarField::Boundary& nuBf = tnu.ref().boundaryFieldRef();
    const volScalarField::Boundary& sigmaBf = strainRate.boundaryField();

    forAll(nuBf, patchi)
    {
        nuBf[patchi] = strainRateFunction_->value(sigmaBf[patchi]);
    }

    return tnu;
}


// ************************************************************************* //
