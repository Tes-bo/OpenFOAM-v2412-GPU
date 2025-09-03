/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
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

#include "nutkWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::nutkWallFunctionFvPatchScalarField::
calcNut() const
{
    const label patchi = patch().index();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
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

    const labelUList& faceCells = patch().faceCells();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    // Viscous sublayer contribution
    const tmp<scalarField> tnutVis = turbModel.nu(patchi);
    const scalarField& nutVis = tnutVis();

    // Calculate y-plus
    const auto yPlus = [&](const label facei) -> scalar
    {
        return (Cmu25*y[facei]*sqrt(k[faceCells[facei]])/nutVis[facei]);
    };

    // Inertial sublayer contribution
    const auto nutLog = [&](const label facei) -> scalar
    {
        const scalar yPlusFace = yPlus(facei);
        return
        (
            nutVis[facei]*yPlusFace*kappa
          / log(max(E*yPlusFace, 1 + 1e-4))
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
                if (yPlus(facei) > yPlusLam)
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
                const scalar yPlusFace = yPlus(facei);
                const scalar Gamma = 0.01*pow4(yPlusFace)/(1 + 5*yPlusFace);
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
                const scalar phiTanh = tanh(pow4(0.1*yPlus(facei)));

                nutw[facei] = phiTanh*b1 + (1 - phiTanh)*b2;
            }
            break;
        }
    }

    nutw -= nutVis;

    return tnutw;
}


void Foam::nutkWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    wallFunctionBlenders()
{}


Foam::nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf)
{}


Foam::nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(4))
{}


Foam::nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    wallFunctionBlenders(wfpsf)
{}


Foam::nutkWallFunctionFvPatchScalarField::nutkWallFunctionFvPatchScalarField
(
    const nutkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    wallFunctionBlenders(wfpsf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::nutkWallFunctionFvPatchScalarField::
yPlus() const
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

    tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();
    tmp<scalarField> tkwc = k.boundaryField()[patchi].patchInternalField();
    const scalarField& kwc = tkwc();

    tmp<scalarField> tnutVis = turbModel.nu(patchi);
    const scalarField& nutVis = tnutVis();

    tmp<scalarField> tnuEff = turbModel.nuEff(patchi);
    const scalarField& nuEff = tnuEff();

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    auto tyPlus = tmp<scalarField>::New(patch().size(), Zero);
    auto& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        // inertial sublayer
        yPlus[facei] = Cmu25*y[facei]*sqrt(kwc[facei])/nutVis[facei];

        if (yPlusLam > yPlus[facei])
        {
            // viscous sublayer
            yPlus[facei] =
                y[facei]*sqrt(nuEff[facei]*magGradUw[facei])/nutVis[facei];
        }
    }

    return tyPlus;
}


void Foam::nutkWallFunctionFvPatchScalarField::write
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
        nutkWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
