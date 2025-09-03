/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2017-2024 OpenCFD Ltd.
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

#include "epsilonWallFunctionFvPatchScalarField.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "turbulenceModel.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::scalar Foam::epsilonWallFunctionFvPatchScalarField::tolerance_ = 1e-5;

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::epsilonWallFunctionFvPatchScalarField::setMaster()
{
    if (master_ != -1)
    {
        return;
    }

    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    label master = -1;
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            if (master == -1)
            {
                master = patchi;
            }

            epf.master() = master;
        }
    }
}


void Foam::epsilonWallFunctionFvPatchScalarField::createAveragingWeights()
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const fvMesh& mesh = epsilon.mesh();

    if (initialised_ && !mesh.changing())
    {
        return;
    }

    volScalarField weights
    (
        IOobject
        (
            "weights",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, Zero)
    );

    DynamicList<label> epsilonPatches(bf.size());
    forAll(bf, patchi)
    {
        if (isA<epsilonWallFunctionFvPatchScalarField>(bf[patchi]))
        {
            epsilonPatches.append(patchi);

            const labelUList& faceCells = bf[patchi].patch().faceCells();
            for (const auto& faceCell : faceCells)
            {
                ++weights[faceCell];
            }
        }
    }

    cornerWeights_.setSize(bf.size());

    for (const auto& patchi : epsilonPatches)
    {
        const fvPatchScalarField& wf = weights.boundaryField()[patchi];
        cornerWeights_[patchi] = 1.0/wf.patchInternalField();
    }

    G_.setSize(internalField().size(), Zero);
    epsilon_.setSize(internalField().size(), Zero);

    initialised_ = true;
}


Foam::epsilonWallFunctionFvPatchScalarField&
Foam::epsilonWallFunctionFvPatchScalarField::epsilonPatch
(
    const label patchi
)
{
    const auto& epsilon =
        static_cast<const volScalarField&>(this->internalField());

    const volScalarField::Boundary& bf = epsilon.boundaryField();

    const auto& epf =
        refCast<const epsilonWallFunctionFvPatchScalarField>(bf[patchi]);

    return const_cast<epsilonWallFunctionFvPatchScalarField&>(epf);
}


void Foam::epsilonWallFunctionFvPatchScalarField::calculateTurbulenceFields
(
    const turbulenceModel& turbulence,
    scalarField& G0,
    scalarField& epsilon0
)
{
    // Accumulate all of the G and epsilon contributions
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            const List<scalar>& w = cornerWeights_[patchi];

            epf.calculate(turbulence, w, epf.patch(), G0, epsilon0);
        }
    }

    // Apply zero-gradient condition for epsilon
    forAll(cornerWeights_, patchi)
    {
        if (!cornerWeights_[patchi].empty())
        {
            epsilonWallFunctionFvPatchScalarField& epf = epsilonPatch(patchi);

            epf == scalarField(epsilon0, epf.patch().faceCells());
        }
    }
}


void Foam::epsilonWallFunctionFvPatchScalarField::calculate
(
    const turbulenceModel& turbModel,
    const List<scalar>& cornerWeights,
    const fvPatch& patch,
    scalarField& G0,
    scalarField& epsilon0
)
{
    const label patchi = patch.index();

    const scalar Cmu25 = pow025(wallCoeffs_.Cmu());
    const scalar Cmu75 = pow(wallCoeffs_.Cmu(), 0.75);
    const scalar kappa = wallCoeffs_.kappa();
    const scalar yPlusLam = wallCoeffs_.yPlusLam();

    const scalarField& y = turbModel.y()[patchi];

    const labelUList& faceCells = patch.faceCells();

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();

    // Calculate y-plus
    const auto yPlus = [&](const label facei) -> scalar
    {
        return
        (
            Cmu25*y[facei]*sqrt(k[faceCells[facei]])/nuw[facei]
        );
    };

    // Contribution from the viscous sublayer
    const auto epsilonVis = [&](const label facei) -> scalar
    {
        return
        (
            2.0*k[faceCells[facei]]*nuw[facei]
          / sqr(y[facei])
        );
    };

    // Contribution from the inertial sublayer
    const auto epsilonLog = [&](const label facei) -> scalar
    {
        return
        (
            Cmu75*pow(k[faceCells[facei]], 1.5)
          / (kappa*y[facei])
        );
    };

    switch (blender_)
    {
        case blenderType::STEPWISE:
        {
            forAll(faceCells, facei)
            {
                if (lowReCorrection_ && yPlus(facei) < yPlusLam)
                {
                    epsilon0[faceCells[facei]] +=
                        cornerWeights[facei]
                      * epsilonVis(facei);
                }
                else
                {
                    epsilon0[faceCells[facei]] +=
                        cornerWeights[facei]
                      * epsilonLog(facei);
                }
            }
            break;
        }

        case blenderType::BINOMIAL:
        {
            forAll(faceCells, facei)
            {
                // (ME:Eqs. 15-16)
                epsilon0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * pow
                    (
                        pow(epsilonVis(facei), n_) + pow(epsilonLog(facei), n_),
                        scalar(1)/n_
                    );
            }
            break;
        }

        case blenderType::MAX:
        {
            forAll(faceCells, facei)
            {
                // (PH:Eq. 27)
                epsilon0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * max(epsilonVis(facei), epsilonLog(facei));
            }
            break;
        }

        case blenderType::EXPONENTIAL:
        {
            forAll(faceCells, facei)
            {
                // (PH:p. 193)
                const scalar yPlusFace = yPlus(facei);
                const scalar Gamma =
                    0.001*pow4(yPlusFace)/(scalar(1) + yPlusFace);
                const scalar invGamma = scalar(1)/(Gamma + ROOTVSMALL);

                epsilon0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * (
                        epsilonVis(facei)*exp(-Gamma)
                      + epsilonLog(facei)*exp(-invGamma)
                    );
            }
            break;
        }

        case blenderType::TANH:
        {
            forAll(faceCells, facei)
            {
                // (KAS:Eqs. 33-34)
                const scalar epsilonVisFace = epsilonVis(facei);
                const scalar epsilonLogFace = epsilonLog(facei);
                const scalar b1 = epsilonVisFace + epsilonLogFace;
                const scalar b2 =
                    pow
                    (
                        pow(epsilonVisFace, 1.2) + pow(epsilonLogFace, 1.2),
                        1.0/1.2
                    );
                const scalar phiTanh = tanh(pow4(0.1*yPlus(facei)));

                epsilon0[faceCells[facei]] +=
                    cornerWeights[facei]
                  * (phiTanh*b1 + (1 - phiTanh)*b2);
            }
            break;
        }
    }

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magGradUw(mag(Uw.snGrad()));

    const tmp<scalarField> tnutw = turbModel.nut(patchi);
    const scalarField& nutw = tnutw();

    forAll(faceCells, facei)
    {
        if (!lowReCorrection_ || (yPlus(facei) > yPlusLam))
        {
            G0[faceCells[facei]] +=
                cornerWeights[facei]
               *(nutw[facei] + nuw[facei])
               *magGradUw[facei]
               *Cmu25*sqrt(k[faceCells[facei]])
               /(kappa*y[facei]);
        }
    }
}


void Foam::epsilonWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    wallFunctionBlenders::writeEntries(os);
    os.writeEntryIfDifferent<bool>("lowReCorrection", false, lowReCorrection_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    wallFunctionBlenders(),
    lowReCorrection_(false),
    initialised_(false),
    master_(-1),
    wallCoeffs_(),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    wallFunctionBlenders(ptf),
    lowReCorrection_(ptf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ptf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
    wallFunctionBlenders(dict, blenderType::STEPWISE, scalar(2)),
    lowReCorrection_(dict.getOrDefault("lowReCorrection", false)),
    initialised_(false),
    master_(-1),
    wallCoeffs_(dict),
    G_(),
    epsilon_(),
    cornerWeights_()
{
    // Apply zero-gradient condition on start-up
    this->extrapolateInternal();
}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf
)
:
    fixedValueFvPatchField<scalar>(ewfpsf),
    wallFunctionBlenders(ewfpsf),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


Foam::epsilonWallFunctionFvPatchScalarField::
epsilonWallFunctionFvPatchScalarField
(
    const epsilonWallFunctionFvPatchScalarField& ewfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ewfpsf, iF),
    wallFunctionBlenders(ewfpsf),
    lowReCorrection_(ewfpsf.lowReCorrection_),
    initialised_(false),
    master_(-1),
    wallCoeffs_(ewfpsf.wallCoeffs_),
    G_(),
    epsilon_(),
    cornerWeights_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalarField& Foam::epsilonWallFunctionFvPatchScalarField::G
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            G_ = 0.0;
        }

        return G_;
    }

    return epsilonPatch(master_).G();
}


Foam::scalarField& Foam::epsilonWallFunctionFvPatchScalarField::epsilon
(
    bool init
)
{
    if (patch().index() == master_)
    {
        if (init)
        {
            epsilon_ = 0.0;
        }

        return epsilon_;
    }

    return epsilonPatch(master_).epsilon(init);
}


void Foam::epsilonWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    forAll(*this, facei)
    {
        const label celli = patch().faceCells()[facei];

        G[celli] = G0[celli];
        epsilon[celli] = epsilon0[celli];
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonWallFunctionFvPatchScalarField::updateWeightedCoeffs
(
    const scalarField& weights
)
{
    if (updated())
    {
        return;
    }

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    setMaster();

    if (patch().index() == master_)
    {
        createAveragingWeights();
        calculateTurbulenceFields(turbModel, G(true), epsilon(true));
    }

    const scalarField& G0 = this->G();
    const scalarField& epsilon0 = this->epsilon();

    typedef DimensionedField<scalar, volMesh> FieldType;

    FieldType& G = db().lookupObjectRef<FieldType>(turbModel.GName());

    FieldType& epsilon = const_cast<FieldType&>(internalField());

    scalarField& epsilonf = *this;

    // Only set the values if the weights are > tolerance
    forAll(weights, facei)
    {
        const scalar w = weights[facei];

        if (w > tolerance_)
        {
            const label celli = patch().faceCells()[facei];

            G[celli] = (1.0 - w)*G[celli] + w*G0[celli];
            epsilon[celli] = (1.0 - w)*epsilon[celli] + w*epsilon0[celli];
            epsilonf[facei] = epsilon[celli];
        }
    }

    fvPatchField<scalar>::updateCoeffs();
}


void Foam::epsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    matrix.setValues(patch().faceCells(), patchInternalField());

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix,
    const Field<scalar>& weights
)
{
    if (manipulatedMatrix())
    {
        return;
    }

    DynamicList<label> constraintCells(weights.size());
    DynamicList<scalar> constraintValues(weights.size());
    const labelUList& faceCells = patch().faceCells();

    const DimensionedField<scalar, volMesh>& fld = internalField();

    forAll(weights, facei)
    {
        // Only set the values if the weights are > tolerance
        if (weights[facei] > tolerance_)
        {
            const label celli = faceCells[facei];

            constraintCells.append(celli);
            constraintValues.append(fld[celli]);
        }
    }

    if (debug)
    {
        Pout<< "Patch: " << patch().name()
            << ": number of constrained cells = " << constraintCells.size()
            << " out of " << patch().size()
            << endl;
    }

    matrix.setValues(constraintCells, constraintValues);

    fvPatchField<scalar>::manipulateMatrix(matrix);
}


void Foam::epsilonWallFunctionFvPatchScalarField::write
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
        epsilonWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
