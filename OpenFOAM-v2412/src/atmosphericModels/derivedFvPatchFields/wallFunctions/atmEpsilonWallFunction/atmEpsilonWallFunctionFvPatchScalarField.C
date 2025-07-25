/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 ENERCON GmbH
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "atmEpsilonWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::atmEpsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    const scalarField& y = turbModel.y()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

    const scalarField magGradUw(mag(Uw.snGrad()));

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar Cmu75 = pow(wallCoeffs_.Cmu(), 0.75);
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    const scalar t = db().time().timeOutputValue();
    const scalarField z0(z0_->value(t));

    #ifdef FULLDEBUG
    for (const auto& z : z0)
    {
        if (z < VSMALL)
        {
            FatalErrorInFunction
                << "z0 field can only contain positive values. "
                << "Please check input field z0."
                << exit(FatalError);
        }
    }
    #endif

    const labelUList& faceCells = patch.faceCells();

    // Set epsilon and G
    forAll(nutw, facei)
    {
        const label celli = faceCells[facei];

        const scalar yPlus = Cmu25*y[facei]*sqrt(k[celli])/nuw[facei];

        const scalar w = cornerWeights[facei];

        // (PGVB:Eq. 7, RH:Eq. 8)
        scalar epsilonc =
            w*Cmu75*pow(k[celli], 1.5)/(kappa*(y[facei] + z0[facei]));

        scalar Gc =
            w
           *(nutw[facei] + nuw[facei])
           *magGradUw[facei]
           *Cmu25*sqrt(k[celli])
           /(kappa*(y[facei] + z0[facei]));

        if (lowReCorrection_ && yPlus < yPlusLam)
        {
            epsilonc = w*2.0*k[celli]*nuw[facei]/sqr(y[facei] + z0[facei]);
            Gc = 0;
        }

        epsilon0[celli] += epsilonc;

        G0[celli] += Gc;
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<bool>("lowReCorrection", false, lowReCorrection_);

    if (z0_)
    {
        z0_->writeData(os);
    }

    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF),
    z0_(nullptr)
{}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    epsilonWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    z0_(ptf.z0_.clone(p.patch()))
{}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    epsilonWallFunctionFvPatchScalarField(p, iF, dict),
    z0_(PatchFunction1<scalar>::New(p.patch(), "z0", dict))
{}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf),
    z0_(ewfpsf.z0_.clone(this->patch().patch()))
{}


Foam::atmEpsilonWallFunctionFvPatchScalarField::
atmEpsilonWallFunctionFvPatchScalarField
(
    const atmEpsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    epsilonWallFunctionFvPatchScalarField(ewfpsf, iF),
    z0_(ewfpsf.z0_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::atmEpsilonWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    epsilonWallFunctionFvPatchScalarField::autoMap(m);

    if (z0_)
    {
        z0_->autoMap(m);
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    epsilonWallFunctionFvPatchScalarField::rmap(ptf, addr);

    const auto& atmpsf =
        refCast<const atmEpsilonWallFunctionFvPatchScalarField>(ptf);
    if (z0_)
    {
        z0_->rmap(atmpsf.z0_(), addr);
    }
}


void Foam::atmEpsilonWallFunctionFvPatchScalarField::write
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
        atmEpsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
