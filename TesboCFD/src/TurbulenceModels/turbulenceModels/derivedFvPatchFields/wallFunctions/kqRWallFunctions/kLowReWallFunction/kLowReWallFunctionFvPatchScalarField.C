/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "kLowReWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::kLowReWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<scalar>("Ceps2", 1.9, Ceps2_);
    os.writeEntryIfDifferent<scalar>("Ck", -0.416, Ck_);
    os.writeEntryIfDifferent<scalar>("Bk", 8.366, Bk_);
    os.writeEntryIfDifferent<scalar>("C", 11.0, C_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kLowReWallFunctionFvPatchScalarField::kLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    Ceps2_(1.9),
    Ck_(-0.416),
    Bk_(8.366),
    C_(11.0),
    wallCoeffs_()
{}


Foam::kLowReWallFunctionFvPatchScalarField::kLowReWallFunctionFvPatchScalarField
(
    const kLowReWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    Ceps2_(ptf.Ceps2_),
    Ck_(ptf.Ck_),
    Bk_(ptf.Bk_),
    C_(ptf.C_),
    wallCoeffs_(ptf.wallCoeffs_)
{}


Foam::kLowReWallFunctionFvPatchScalarField::kLowReWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    Ceps2_
    (
        dict.getCheckOrDefault<scalar>
        (
            "Ceps2",
            1.9,
            scalarMinMax::ge(SMALL)
        )
    ),
    Ck_(dict.getOrDefault<scalar>("Ck", -0.416)),
    Bk_(dict.getOrDefault<scalar>("Bk", 8.366)),
    C_(dict.getOrDefault<scalar>("C", 11.0)),
    wallCoeffs_(dict)
{}


Foam::kLowReWallFunctionFvPatchScalarField::kLowReWallFunctionFvPatchScalarField
(
    const kLowReWallFunctionFvPatchScalarField& kwfpsf
)
:
    fixedValueFvPatchField<scalar>(kwfpsf),
    Ceps2_(kwfpsf.Ceps2_),
    Ck_(kwfpsf.Ck_),
    Bk_(kwfpsf.Bk_),
    C_(kwfpsf.C_),
    wallCoeffs_(kwfpsf.wallCoeffs_)
{}


Foam::kLowReWallFunctionFvPatchScalarField::kLowReWallFunctionFvPatchScalarField
(
    const kLowReWallFunctionFvPatchScalarField& kwfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(kwfpsf, iF),
    Ceps2_(kwfpsf.Ceps2_),
    Ck_(kwfpsf.Ck_),
    Bk_(kwfpsf.Bk_),
    C_(kwfpsf.C_),
    wallCoeffs_(kwfpsf.wallCoeffs_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kLowReWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

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

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    scalarField& kw = *this;

    // Set k wall values
    forAll(kw, facei)
    {
        const label celli = patch().faceCells()[facei];
        const scalar uTau = Cmu25*sqrt(k[celli]);
        const scalar yPlus = uTau*y[facei]/nuw[facei];

        if (yPlus > yPlusLam)
        {
            kw[facei] = Ck_/kappa*log(yPlus) + Bk_;
        }
        else
        {
            const scalar Cf =
                1.0/sqr(yPlus + C_) + 2.0*yPlus/pow3(C_) - 1.0/sqr(C_);
            kw[facei] = 2400.0/sqr(Ceps2_)*Cf;
        }

        kw[facei] *= sqr(uTau);
    }

    // Limit kw to avoid failure of the turbulence model due to division by kw
    kw = max(kw, SMALL);

    fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void Foam::kLowReWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        kLowReWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
