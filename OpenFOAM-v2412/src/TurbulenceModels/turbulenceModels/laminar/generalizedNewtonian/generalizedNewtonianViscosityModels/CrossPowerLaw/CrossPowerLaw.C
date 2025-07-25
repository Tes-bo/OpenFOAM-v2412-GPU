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

#include "CrossPowerLaw.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarModels
{
namespace generalizedNewtonianViscosityModels
{
    defineTypeNameAndDebug(CrossPowerLaw, 0);

    addToRunTimeSelectionTable
    (
        generalizedNewtonianViscosityModel,
        CrossPowerLaw,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
CrossPowerLaw
(
    const dictionary& viscosityProperties
)
:
    generalizedNewtonianViscosityModel(viscosityProperties),
    CrossPowerLawCoeffs_
    (
        viscosityProperties.optionalSubDict(typeName + "Coeffs")
    ),
    nuInf_("nuInf", dimViscosity, CrossPowerLawCoeffs_),
    m_("m", dimTime, CrossPowerLawCoeffs_),
    n_("n", dimless, CrossPowerLawCoeffs_),
    tauStar_("tauStar", dimViscosity/dimTime, CrossPowerLawCoeffs_)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
read
(
    const dictionary& viscosityProperties
)
{
    generalizedNewtonianViscosityModel::read(viscosityProperties);

    CrossPowerLawCoeffs_ =
        viscosityProperties.optionalSubDict(typeName + "Coeffs");

    CrossPowerLawCoeffs_.readEntry("nuInf", nuInf_);
    CrossPowerLawCoeffs_.readEntry("m", m_);
    CrossPowerLawCoeffs_.readEntry("n", n_);

    return true;
}


Foam::tmp<Foam::volScalarField>
Foam::laminarModels::generalizedNewtonianViscosityModels::CrossPowerLaw::
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
      /
        (
            scalar(1)
          + pow
            (
                (tauStar_.value() > 0)
              ? nu0*strainRate/tauStar_
              : m_*strainRate,
                n_
            )
        )
    );
}


// ************************************************************************* //
