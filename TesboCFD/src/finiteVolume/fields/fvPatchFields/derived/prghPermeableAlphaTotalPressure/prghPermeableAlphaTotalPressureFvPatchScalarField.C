/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "prghPermeableAlphaTotalPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "gravityMeshObject.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::
prghPermeableAlphaTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(p, iF),
    p0_(nullptr),
    phiName_("phi"),
    rhoName_("rho"),
    UName_("U"),
    alphaName_("none"),
    alphaMin_(1.0),
    curTimeIndex_(-1)
{
    refValue() = 0.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;
}


Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::
prghPermeableAlphaTotalPressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<scalar>(p, iF),
    p0_(PatchFunction1<scalar>::New(p.patch(), "p", dict)),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    UName_(dict.getOrDefault<word>("U", "U")),
    alphaName_(dict.getOrDefault<word>("alpha", "none")),
    alphaMin_(dict.getOrDefault<scalar>("alphaMin", 1)),
    curTimeIndex_(-1)
{
    refValue() = 1.0;
    refGrad() = 0.0;
    valueFraction() = 0.0;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<scalar>::operator=(refValue());
    }
}


Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::
prghPermeableAlphaTotalPressureFvPatchScalarField
(
    const prghPermeableAlphaTotalPressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<scalar>(ptf, p, iF, mapper),
    p0_(ptf.p0_.clone(p.patch())),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_),
    UName_(ptf.UName_),
    alphaName_(ptf.alphaName_),
    alphaMin_(ptf.alphaMin_),
    curTimeIndex_(-1)
{}


Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::
prghPermeableAlphaTotalPressureFvPatchScalarField
(
    const prghPermeableAlphaTotalPressureFvPatchScalarField& tppsf
)
:
    mixedFvPatchField<scalar>(tppsf),
    p0_(tppsf.p0_.clone(this->patch().patch())),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    UName_(tppsf.UName_),
    alphaName_(tppsf.alphaName_),
    alphaMin_(tppsf.alphaMin_),
    curTimeIndex_(-1)
{}


Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::
prghPermeableAlphaTotalPressureFvPatchScalarField
(
    const prghPermeableAlphaTotalPressureFvPatchScalarField& tppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField<scalar>(tppsf, iF),
    p0_(tppsf.p0_.clone(this->patch().patch())),
    phiName_(tppsf.phiName_),
    rhoName_(tppsf.rhoName_),
    UName_(tppsf.UName_),
    alphaName_(tppsf.alphaName_),
    alphaMin_(tppsf.alphaMin_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchField<scalar>::autoMap(m);

    if (p0_)
    {
        p0_->autoMap(m);
    }
}


void Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchField<scalar>::rmap(ptf, addr);

    const auto& tptf =
        refCast<const prghPermeableAlphaTotalPressureFvPatchScalarField>(ptf);

    if (p0_)
    {
        p0_->rmap(tptf.p0_(), addr);
    }
}


void Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::updateSnGrad
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }

    const scalarField& rhop =
        patch().lookupPatchField<volScalarField>(rhoName_);

    const scalarField& phip =
        patch().lookupPatchField<surfaceScalarField>(phiName_);

    const vectorField& Up =
        patch().lookupPatchField<volVectorField>(UName_);

    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(db().time());

    const auto& hRef =
        db().lookupObject<uniformDimensionedScalarField>("hRef");

    const dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRef
      : dimensionedScalar(g.dimensions()*dimLength, 0)
    );

    const scalar t = db().time().timeOutputValue();

    tmp<scalarField> p
    (
        p0_->value(t)
      - 0.5*rhop*(neg(phip))*magSqr(Up)
      - rhop*((g.value() & patch().Cf()) - ghRef.value())
    );

    refValue() = p;

    refGrad() = snGradp;

    if (alphaName_ != "none")
    {
        const scalarField& alphap =
            patch().lookupPatchField<volScalarField>(alphaName_);

        tmp<scalarField> alphaCut(pos(alphap - alphaMin_));
        valueFraction() = 1 - alphaCut;
    }

    if (debug)
    {
        const scalar phi = gSum(-phip);
        Info<< valueFraction() << endl;
        Info<< patch().boundaryMesh().mesh().name() << ':'
            << patch().name() << ':'
            << this->internalField().name() << " :"
            << " mass flux[Kg/s]:" << phi
            << endl;
    }

    curTimeIndex_ = this->db().time().timeIndex();

    mixedFvPatchField<scalar>::updateCoeffs();
}


void Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }
}


void Foam::prghPermeableAlphaTotalPressureFvPatchScalarField::write
(
    Ostream& os
) const
{
    mixedFvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntryIfDifferent<word>("U", "U", UName_);
    os.writeEntryIfDifferent<word>("alpha", "none", alphaName_);
    os.writeEntryIfDifferent<scalar>("alphaMin", 1, alphaMin_);

    if (p0_)
    {
        p0_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        prghPermeableAlphaTotalPressureFvPatchScalarField
    );

}

// ************************************************************************* //
