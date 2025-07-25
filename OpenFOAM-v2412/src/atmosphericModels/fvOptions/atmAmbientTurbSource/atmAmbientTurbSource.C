/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "atmAmbientTurbSource.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmAmbientTurbSource, 0);
    addToRunTimeSelectionTable(option, atmAmbientTurbSource, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmAmbientTurbSource::atmAmbientTurbSource
(
    const word& sourceName,
    const word& modelType,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    fv::cellSetOption(sourceName, modelType, dict, mesh),
    isEpsilon_(false),
    rhoName_(coeffs_.getOrDefault<word>("rho", "rho")),
    kAmb_
    (
        dimensionedScalar
        (
            sqr(dimLength)/sqr(dimTime),
            coeffs_.getCheckOrDefault<scalar>
            (
                "kAmb",
                SMALL,
                [&](const scalar k){ return k > -VSMALL; }
            )
        )
    ),
    epsilonAmb_
    (
        dimensionedScalar
        (
            sqr(dimLength)/pow3(dimTime),
            coeffs_.getOrDefault<scalar>
            (
                "epsilonAmb",
                Zero
            )
        )
    ),
    omegaAmb_
    (
        dimensionedScalar
        (
            dimless/dimTime,
            coeffs_.getOrDefault<scalar>
            (
                "omegaAmb",
                Zero
            )
        )
    ),
    Cmu_(Zero),
    C2_(Zero)
{
    const auto* turbPtr =
        mesh_.findObject<turbulenceModel>
        (
            turbulenceModel::propertiesName
        );

    if (!turbPtr)
    {
        FatalErrorInFunction
            << "Unable to find a turbulence model."
            << abort(FatalError);
    }

    fieldNames_.resize(2);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (tepsilon.is_reference())
    {
        isEpsilon_ = true;
        fieldNames_[0] = tepsilon().name();

        const dictionary& turbDict = turbPtr->coeffDict();

        C2_.read("C2", turbDict);
    }
    else if (tomega.is_reference())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();

        const dictionary& turbDict = turbPtr->coeffDict();

        Cmu_.read("betaStar", turbDict);
        C2_.read("C2", turbDict);
    }
    else
    {
        FatalErrorInFunction
            << "Needs either epsilon or omega field."
            << abort(FatalError);
    }

    fieldNames_[1] = turbPtr->k()().name();

    fv::option::resetApplied();

    Log << "    Applying atmAmbientTurbSource to: "
        << fieldNames_[0] << " and " << fieldNames_[1]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmAmbientTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmAmbientTurbSourceK
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
        return;
    }

    if (isEpsilon_)
    {
        atmAmbientTurbSourceEpsilon
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
    else
    {
        atmAmbientTurbSourceOmega
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
}


void Foam::fv::atmAmbientTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmAmbientTurbSourceK(geometricOneField(), rho, eqn, fieldi);
        return;
    }

    if (isEpsilon_)
    {
        atmAmbientTurbSourceEpsilon(geometricOneField(), rho, eqn, fieldi);
    }
    else
    {
        atmAmbientTurbSourceOmega(geometricOneField(), rho, eqn, fieldi);
    }
}


void Foam::fv::atmAmbientTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (fieldi == 1)
    {
        atmAmbientTurbSourceK(alpha, rho, eqn, fieldi);
        return;
    }

    if (isEpsilon_)
    {
        atmAmbientTurbSourceEpsilon(alpha, rho, eqn, fieldi);
    }
    else
    {
        atmAmbientTurbSourceOmega(alpha, rho, eqn, fieldi);
    }
}


// ************************************************************************* //
