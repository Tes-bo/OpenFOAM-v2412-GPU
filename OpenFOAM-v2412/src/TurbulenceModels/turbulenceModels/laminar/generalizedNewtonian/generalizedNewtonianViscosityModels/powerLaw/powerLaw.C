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

#include "powerLaw.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(powerLaw, 0);

    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        powerLaw,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::powerLaw::powerLaw
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    powerLawCoeffs_(viscosityProperties.optionalSubDict(typeName + "Coeffs")),
    n_("n", dimless, powerLawCoeffs_),
    nuMin_("nuMin", dimViscosity, powerLawCoeffs_),
    nuMax_("nuMax", dimViscosity, powerLawCoeffs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::powerLaw::read
(
    const dictionary& viscosityProperties
)
{
    generalizedNewtonianViscosityModel::read(viscosityProperties);

    powerLawCoeffs_ = viscosityProperties.optionalSubDict(typeName + "Coeffs");

    powerLawCoeffs_.readEntry("n", n_);
    powerLawCoeffs_.readEntry("nuMin", nuMin_);
    powerLawCoeffs_.readEntry("nuMax", nuMax_);

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::powerLaw::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return max
    (
        nuMin_,
        min
        (
            nuMax_,
            nu0*pow
            (
                max
                (
                    dimensionedScalar("one", dimTime, 1)*strainRate,
                    dimensionedScalar("small", dimless, SMALL)
                ),
                n_.value() - scalar(1)
            )
        )
    );
}


// ************************************************************************* //
