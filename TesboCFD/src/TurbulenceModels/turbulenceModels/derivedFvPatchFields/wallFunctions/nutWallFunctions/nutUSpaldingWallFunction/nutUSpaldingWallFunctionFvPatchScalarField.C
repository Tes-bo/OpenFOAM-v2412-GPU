/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "nutUSpaldingWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUSpaldingWallFunctionFvPatchScalarField::calcNut() const
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

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magGradU(mag(Uw.snGrad()));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();


    // Calculate new viscosity
    tmp<scalarField> tnutw
    (
        max
        (
            scalar(0),
            sqr(calcUTau(magGradU))/(magGradU + ROOTVSMALL) - nuw
        )
    );

    if (tolerance_ != 0.01)
    {
        // User-specified tolerance. Find out if current nut already satisfies
        // eqns.

        // Run ut for one iteration
        scalarField err;
        tmp<scalarField> UTau(calcUTau(magGradU, 1, err));

        // Preserve nutw if the initial error (err) already lower than the
        // tolerance.

        scalarField& nutw = tnutw.ref();
        forAll(err, facei)
        {
            if (err[facei] < tolerance_)
            {
                nutw[facei] = this->operator[](facei);
            }
        }
    }
    return tnutw;
}


Foam::tmp<Foam::scalarField>
Foam::nutUSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU
) const
{
    scalarField err;
    return calcUTau(magGradU, maxIter_, err);
}


Foam::tmp<Foam::scalarField>
Foam::nutUSpaldingWallFunctionFvPatchScalarField::calcUTau
(
    const scalarField& magGradU,
    const label maxIter,
    scalarField& err
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

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const scalar kappa = wallCoeffs_.kappa();
    const scalar E = wallCoeffs_.E();

    const scalarField& nutw = *this;

    auto tuTau = tmp<scalarField>::New(patch().size(), Zero);
    auto& uTau = tuTau.ref();

    err.setSize(uTau.size());
    err = 0.0;

    forAll(uTau, facei)
    {
        scalar ut = sqrt((nutw[facei] + nuw[facei])*magGradU[facei]);
        // Note: for exact restart seed with laminar viscosity only:
        //scalar ut = sqrt(nuw[facei]*magGradU[facei]);

        if (ROOTVSMALL < ut)
        {
            int iter = 0;

            do
            {
                const scalar kUu = min(kappa*magUp[facei]/ut, scalar(50));
                const scalar fkUu = exp(kUu) - 1 - kUu*(1 + 0.5*kUu);

                const scalar f =
                    - ut*y[facei]/nuw[facei]
                    + magUp[facei]/ut
                    + 1.0/E*(fkUu - 1.0/6.0*kUu*sqr(kUu));

                const scalar df =
                    y[facei]/nuw[facei]
                  + magUp[facei]/sqr(ut)
                  + 1.0/E*kUu*fkUu/ut;

                const scalar uTauNew = ut + f/df;
                err[facei] = mag((ut - uTauNew)/ut);
                ut = uTauNew;

                //iterations_++;

            } while
            (
                ut > ROOTVSMALL
             && err[facei] > tolerance_
             && ++iter < maxIter
            );

            uTau[facei] = max(scalar(0), ut);

            //invocations_++;
            //if (iter > 1)
            //{
            //    nontrivial_++;
            //}
            //if (iter >= maxIter_)
            //{
            //    nonconvergence_++;
            //}
        }
    }

    return tuTau;
}


void Foam::nutUSpaldingWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<label>("maxIter", 10, maxIter_);
    os.writeEntryIfDifferent<scalar>("tolerance", 0.01, tolerance_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(p, iF),
    maxIter_(10),
    tolerance_(0.01)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
    maxIter_(ptf.maxIter_),
    tolerance_(ptf.tolerance_)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    maxIter_(dict.getOrDefault<label>("maxIter", 10)),
    tolerance_(dict.getOrDefault<scalar>("tolerance", 0.01))
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


Foam::nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& wfpsf
)
:
    nutWallFunctionFvPatchScalarField(wfpsf),
    maxIter_(wfpsf.maxIter_),
    tolerance_(wfpsf.tolerance_)
    //invocations_(wfpsf.invocations_),
    //nontrivial_(wfpsf.nontrivial_),
    //nonconvergence_(wfpsf.nonconvergence_),
    //iterations_(wfpsf.iterations_)
{}


Foam::nutUSpaldingWallFunctionFvPatchScalarField::
nutUSpaldingWallFunctionFvPatchScalarField
(
    const nutUSpaldingWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(wfpsf, iF),
    maxIter_(wfpsf.maxIter_),
    tolerance_(wfpsf.tolerance_)
    //invocations_(0),
    //nontrivial_(0),
    //nonconvergence_(0),
    //iterations_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::nutUSpaldingWallFunctionFvPatchScalarField::
~nutUSpaldingWallFunctionFvPatchScalarField()
{
    //if (debug)
    //{
    //    Info<< "nutUSpaldingWallFunctionFvPatchScalarField :"
    //        << " total invocations:"
    //        << returnReduce(invocations_, sumOp<label>())
    //        << " total iterations:"
    //        << returnReduce(iterations_, sumOp<label>())
    //        << " total non-convergence:"
    //        << returnReduce(nonconvergence_, sumOp<label>())
    //        << " total non-trivial:"
    //        << returnReduce(nontrivial_, sumOp<label>())
    //        << endl;
    //}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::nutUSpaldingWallFunctionFvPatchScalarField::yPlus() const
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

    const fvPatchVectorField& Uw = U(turbModel).boundaryField()[patchi];

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    return y*calcUTau(mag(Uw.snGrad()))/nuw;
}


void Foam::nutUSpaldingWallFunctionFvPatchScalarField::write
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
        nutUSpaldingWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
