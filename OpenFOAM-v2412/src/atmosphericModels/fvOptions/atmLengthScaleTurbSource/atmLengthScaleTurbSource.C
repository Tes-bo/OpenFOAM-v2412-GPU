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

#include "atmLengthScaleTurbSource.H"
#include "geometricOneField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{
    defineTypeNameAndDebug(atmLengthScaleTurbSource, 0);
    addToRunTimeSelectionTable(option, atmLengthScaleTurbSource, dictionary);
}
}


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField::Internal>Foam::fv::atmLengthScaleTurbSource::
calcC1Star
(
    const volScalarField::Internal& k,
    const volScalarField::Internal& epsilon
) const
{
    // Mixing-length scale estimation (P:Eq. 10.37 & p. 374)
    tmp<volScalarField::Internal> L_(pow(Cmu_, 0.75)*pow(k, 1.5)/epsilon);

    // (AC:Eq. 16) wherein the exponentiation "n_" is not present.
    // "n_" is an ad-hoc implementation.
    return (C2_ - C1_)*pow(L_/Lmax_, n_);
}


Foam::tmp<Foam::volScalarField::Internal> Foam::fv::atmLengthScaleTurbSource::
calcGammaStar
(
    const volScalarField::Internal& k,
    const volScalarField::Internal& omega,
    const volScalarField::Internal& gamma,
    const volScalarField::Internal& beta
) const
{
    // (L:Eq. 3.20)
    tmp<volScalarField::Internal> L_(sqrt(k)/(pow025(Cmu_)*omega));

    // (L:Eq. 3.34)
    return (gamma - beta)*pow(L_/Lmax_, n_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fv::atmLengthScaleTurbSource::atmLengthScaleTurbSource
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
    Lmax_
    (
        dimensionedScalar
        (
            dimLength,
            coeffs_.getCheckOrDefault<scalar>
            (
                "Lmax",
                41.575,
                [&](const scalar Lmax){ return Lmax > SMALL; }
            )
        )
    ),
    n_
    (
        dimensionedScalar
        (
            dimless,
            coeffs_.getCheckOrDefault<scalar>
            (
                "n",
                3.0,
                [&](const scalar n){ return n > SMALL; }
            )
        )
    ),
    Cmu_(Zero),
    C1_(Zero),
    C2_(Zero),
    C3_(Zero)
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

    fieldNames_.resize(1);

    tmp<volScalarField> tepsilon = turbPtr->epsilon();
    tmp<volScalarField> tomega = turbPtr->omega();

    if (tepsilon.is_reference())
    {
        isEpsilon_ = true;
        fieldNames_[0] = tepsilon().name();

        const dictionary& turbDict = turbPtr->coeffDict();
        Cmu_.read("Cmu", turbDict);
        C1_.read("C1", turbDict);
        C2_.read("C2", turbDict);
        C3_.read("C3", turbDict);
    }
    else if (tomega.is_reference())
    {
        isEpsilon_ = false;
        fieldNames_[0] = tomega().name();

        const dictionary& turbDict = turbPtr->coeffDict();
        Cmu_.read("betaStar", turbDict);
    }
    else
    {
        FatalErrorInFunction
            << "Needs either epsilon or omega field."
            << abort(FatalError);
    }

    fv::option::resetApplied();

    Log << "    Applying atmLengthScaleTurbSource to: " << fieldNames_[0]
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fv::atmLengthScaleTurbSource::addSup
(
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmLengthScaleTurbSourceEpsilon
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
    else
    {
        atmLengthScaleTurbSourceOmega
        (
            geometricOneField(),
            geometricOneField(),
            eqn,
            fieldi
        );
    }
}


void Foam::fv::atmLengthScaleTurbSource::addSup
(
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmLengthScaleTurbSourceEpsilon(geometricOneField(), rho, eqn, fieldi);
    }
    else
    {
        atmLengthScaleTurbSourceOmega(geometricOneField(), rho, eqn, fieldi);
    }
}


void Foam::fv::atmLengthScaleTurbSource::addSup
(
    const volScalarField& alpha,
    const volScalarField& rho,
    fvMatrix<scalar>& eqn,
    const label fieldi
)
{
    if (isEpsilon_)
    {
        atmLengthScaleTurbSourceEpsilon(alpha, rho, eqn, fieldi);
    }
    else
    {
        atmLengthScaleTurbSourceOmega(alpha, rho, eqn, fieldi);
    }
}


// ************************************************************************* //
