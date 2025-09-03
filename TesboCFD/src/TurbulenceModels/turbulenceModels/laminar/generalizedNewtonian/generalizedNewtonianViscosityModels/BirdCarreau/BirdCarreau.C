/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd
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

#include "BirdCarreau.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(BirdCarreau, 0);
    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        BirdCarreau,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
BirdCarreau
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    BirdCarreauCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nuInf_("nuInf", dimViscosity, BirdCarreauCoeffs_),
    k_("k", dimTime, BirdCarreauCoeffs_),
    n_("n", dimless, BirdCarreauCoeffs_),
    a_
    (
        BirdCarreauCoeffs_.getOrDefault
        (
            "a",
            dimensionedScalar("a", dimless, 2)
        )
    )
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
read
(
    const dictionary& viscosityProperties
)
{
    generalizedNewtonianViscosityModel::read(viscosityProperties);

    BirdCarreauCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    BirdCarreauCoeffs_.readEntry("nuInf", nuInf_);
    BirdCarreauCoeffs_.readEntry("k", k_);
    BirdCarreauCoeffs_.readEntry("n", n_);

    a_ = BirdCarreauCoeffs_.getOrDefault
    (
        "a",
        dimensionedScalar("a", dimless, 2)
    );

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::BirdCarreau::
nu
(
    const volScalarField& nu0,
    const volScalarField& strainRate
) const
{
    return
    (
        nuInf_
      + (nu0 - nuInf_)
      * pow
        (
            scalar(1)
          + pow
            (
                k_*strainRate,
                a_
            ),
            (n_ - scalar(1))/a_
        )
    );
}


// ************************************************************************* //
