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

#include "atmAlphatkWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar atmAlphatkWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label atmAlphatkWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void atmAlphatkWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


void atmAlphatkWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<scalar>("Cmu", 0.09, Cmu_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);

    if (Pr_)
    {
        Pr_->writeData(os);
    }
    if (Prt_)
    {
        Prt_->writeData(os);
    }
    if (z0_)
    {
        z0_->writeData(os);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmAlphatkWallFunctionFvPatchScalarField::
atmAlphatkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Cmu_(0.09),
    kappa_(0.41),
    Pr_(nullptr),
    Prt_(nullptr),
    z0_(nullptr)
{
    checkType();
}


atmAlphatkWallFunctionFvPatchScalarField::
atmAlphatkWallFunctionFvPatchScalarField
(
    const atmAlphatkWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Cmu_(ptf.Cmu_),
    kappa_(ptf.kappa_),
    Pr_(ptf.Pr_.clone()),
    Prt_(ptf.Prt_.clone(p.patch())),
    z0_(ptf.z0_.clone(p.patch()))
{
    checkType();
}


atmAlphatkWallFunctionFvPatchScalarField::
atmAlphatkWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Cmu_
    (
        dict.getCheckOrDefault<scalar>
        (
            "Cmu",
            0.09,
            scalarMinMax::ge(SMALL)
        )
    ),
    kappa_
    (
        dict.getCheckOrDefault<scalar>
        (
            "kappa",
            0.41,
            scalarMinMax::ge(SMALL)
        )
    ),
    Pr_(Function1<scalar>::New("Pr", dict, &db())),
    Prt_(PatchFunction1<scalar>::New(p.patch(), "Prt", dict)),
    z0_(PatchFunction1<scalar>::New(p.patch(), "z0", dict))
{
    checkType();
}


atmAlphatkWallFunctionFvPatchScalarField::
atmAlphatkWallFunctionFvPatchScalarField
(
    const atmAlphatkWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    Pr_(wfpsf.Pr_.clone()),
    Prt_(wfpsf.Prt_.clone(this->patch().patch())),
    z0_(wfpsf.z0_.clone(this->patch().patch()))
{
    checkType();
}


atmAlphatkWallFunctionFvPatchScalarField::
atmAlphatkWallFunctionFvPatchScalarField
(
    const atmAlphatkWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Cmu_(wfpsf.Cmu_),
    kappa_(wfpsf.kappa_),
    Pr_(wfpsf.Pr_.clone()),
    Prt_(wfpsf.Prt_.clone(this->patch().patch())),
    z0_(wfpsf.z0_.clone(this->patch().patch()))
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmAlphatkWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model
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

    const scalar Cmu25 = pow025(Cmu_);

    const scalar t = db().time().timeOutputValue();
    const scalar Pr = Pr_->value(t);

    #ifdef FULLDEBUG
    if (Pr < VSMALL)
    {
        FatalErrorInFunction
            << "Pr cannot be negative or zero. "
            << "Please check input Pr = " << Pr
            << exit(FatalError);
    }
    #endif

    const scalarField Prt(Prt_->value(t));
    const scalarField z0(z0_->value(t));

    #ifdef FULLDEBUG
    forAll(Prt, i)
    {
        if (Prt[i] < VSMALL || z0[i] < VSMALL)
        {
            FatalErrorInFunction
                << "Elements of input surface fields can only be positive. "
                << "Please check input fields z0 and Prt."
                << exit(FatalError);
        }
    }
    #endif

    const labelUList& faceCells = patch().faceCells();

    scalarField& alphatw = *this;

    forAll(alphatw, facei)
    {
        const label celli = faceCells[facei];

        const scalar uStar = Cmu25*Foam::sqrt(k[celli]);
        const scalar Edash = (y[facei] + z0[facei])/(z0[facei] + 1e-4);

        // Update turbulent thermal conductivity
        alphatw[facei] =
            uStar*kappa_*y[facei]/(Prt[facei]*log(max(Edash, 1 + 1e-4)))
          + nuw[facei]/Pr;
    }

    // lower bound values to avoid unrealistic
    // negative temperatures on the ground
    alphatw = max(alphatw, scalar(0.01));

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void atmAlphatkWallFunctionFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);

    if (Prt_)
    {
        Prt_->autoMap(m);
    }
    if (z0_)
    {
        z0_->autoMap(m);
    }
}


void atmAlphatkWallFunctionFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);

    const auto& nrwfpsf =
        refCast<const atmAlphatkWallFunctionFvPatchScalarField>(ptf);

    if (Prt_)
    {
        Prt_->rmap(nrwfpsf.Prt_(), addr);
    }
    if (z0_)
    {
        z0_->rmap(nrwfpsf.z0_(), addr);
    }
}


void atmAlphatkWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmAlphatkWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
