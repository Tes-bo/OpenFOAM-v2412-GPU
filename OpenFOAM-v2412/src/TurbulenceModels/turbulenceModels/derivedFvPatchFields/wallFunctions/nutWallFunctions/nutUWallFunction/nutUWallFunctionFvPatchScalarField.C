/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "nutUWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUWallFunctionFvPatchScalarField::calcNut() const
{
    const label patchi = patch().index();

    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    const scalarField& yPlus = tyPlus();

    // Viscous sublayer contribution
    const tmp<scalarField> tnutVis = turbModel.nu(patchi);
    const scalarField& nutVis = tnutVis();

    // Inertial sublayer contribution
    const auto nutLog = [&](const label facei) -> scalar
    {
        const scalar yPlusFace = yPlus[facei];
        return
        (
            nutVis[facei]*yPlusFace*kappa/log(max(E*yPlusFace, 1 + 1e-4))
        );
    };

    auto tnutw = tmp<scalarField>::New(patch().size(), Zero);
    auto& nutw = tnutw.ref();

    switch (blender_)
    {
        case blenderType::STEPWISE:
        {
            forAll(nutw, facei)
            {
                if (yPlus[facei] > yPlusLam)
                {
                    nutw[facei] = nutLog(facei);
                }
                else
                {
                    nutw[facei] = nutVis[facei];
                }
            }
            break;
        }

        case blenderType::MAX:
        {
            forAll(nutw, facei)
            {
                // (PH:Eq. 27)
                nutw[facei] = max(nutVis[facei], nutLog(facei));
            }
            break;
        }

        case blenderType::BINOMIAL:
        {
            forAll(nutw, facei)
            {
                // (ME:Eqs. 15-16)
                nutw[facei] =
                    pow
                    (
                        pow(nutVis[facei], n_) + pow(nutLog(facei), n_),
                        scalar(1)/n_
                    );
            }
            break;
        }

        case blenderType::EXPONENTIAL:
        {
            forAll(nutw, facei)
            {
                // (PH:Eq. 31)
                const scalar Gamma =
                    0.01*pow4(yPlus[facei])/(1 + 5*yPlus[facei]);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                nutw[facei] =
                    nutVis[facei]*exp(-Gamma) + nutLog(facei)*exp(-invGamma);
            }
            break;
        }

        case blenderType::TANH:
        {
            forAll(nutw, facei)
            {
                // (KAS:Eqs. 33-34)
                const scalar nutLogFace = nutLog(facei);
                const scalar b1 = nutVis[facei] + nutLogFace;
                const scalar b2 =
                    pow
                    (
                        pow(nutVis[facei], 1.2) + pow(nutLogFace, 1.2),
                        1.0/1.2
                    );
                const scalar phiTanh = tanh(pow4(0.1*yPlus[facei]));

                nutw[facei] = phiTanh*b1 + (1 - phiTanh)*b2;
            }
            break;
        }
    }

    nutw -= nutVis;

    return tnutw;
}


Foam::tmp<Foam::scalarField>
Foam::nutUWallFunctionFvPatchScalarField::calcYPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    auto tyPlus = tmp<scalarField>::New(patch().size(), Zero);
    auto& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        const scalar kappaRe = kappa*magUp[facei]*y[facei]/nuw[facei];

        scalar yp = yPlusLam;
        const scalar ryPlusLam = 1.0/yp;

        int iter = 0;
        scalar yPlusLast = 0.0;

        do
        {
            yPlusLast = yp;
            yp = (kappaRe + yp)/(1.0 + log(E*yp));

        } while (mag(ryPlusLam*(yp - yPlusLast)) > 0.01 && ++iter < 10 );

        yPlus[facei] = max(scalar(0), yp);
    }

    return tyPlus;
}


void Foam::nutUWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    wallFunctionBlenders()
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(4))
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf),
    wallFunctionBlenders(sawfpsf)
{}


Foam::nutUWallFunctionFvPatchScalarField::nutUWallFunctionFvPatchScalarField
(
    const nutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf, iF),
    wallFunctionBlenders(sawfpsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField& y = turbModel.y()[patchi];

    tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    tmp<scalarField> tnuEff = turbModel.nuEff(patchi);
    const scalarField& nuEff = tnuEff();

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    tmp<scalarField> tyPlus = calcYPlus(magUp);
    scalarField& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        if (yPlusLam > yPlus[facei])
        {
            // viscous sublayer
            yPlus[facei] =
                y[facei]*sqrt(nuEff[facei]*magGradUw[facei])/nuw[facei];
        }
    }

    return tyPlus;
}


void Foam::nutUWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    nutWallFunctionFvPatchScalarField::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        nutUWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
