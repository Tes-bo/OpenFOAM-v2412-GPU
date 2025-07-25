/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "turbulentDigitalFilterInletFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "faceAreaWeightAMI.H"
#include "turbulentDFSEMInletFvPatchVectorField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::vector
Foam::turbulentDigitalFilterInletFvPatchField<Type>::calcPatchNormal() const
{
    const vectorField nf(this->patch().nf());

    // Patch normal points into domain
    vector patchNormal(-gAverage(nf));

    // Check that patch is planar
    const scalar error = max(magSqr(patchNormal + nf));

    if (error > SMALL)
    {
        WarningInFunction
            << "Patch " << this->patch().name() << " is not planar"
            << endl;
    }

    return patchNormal.normalise();
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::initialisePatch()
{
    L_.initialise();

    AMIPtr_->calculate(this->patch().patch(), L_.patch());

    patchNormal_ = calcPatchNormal();
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::mapL
(
    Field<Type>& fld
)
{
    Field<Type> sourceFld;

    if (Pstream::master())
    {
        sourceFld = L_.convolve();
        L_.shift();
        L_.refill();
    }

    // Map two-point correlations (integral scales)
    plusEqOp<Type> cop;
    AMIPtr_->interpolateToSource
    (
        sourceFld,
        multiplyWeightedOp<Type, plusEqOp<Type>>(cop),
        fld,
        UList<Type>::null()
    );

    // Map forward-stepwise method correlations if requested
    if (L_.fsm())
    {
        L_.correlate(fld);
    }
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::mapR
(
    scalarField& fld
) const
{
    const scalar t = this->db().time().timeOutputValue();
    scalarField R(Rptr_->value(t));

    // Lund-Wu-Squires transformation for scalar fields
    R = Foam::sqrt(R);

    // Map transformed Reynolds stresses field onto patch for scalar fields
    fld *= R;
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::mapR
(
    vectorField& fld
) const
{
    const scalar t = this->db().time().timeOutputValue();
    symmTensorField R(Rptr_->value(t));

    // Lund-Wu-Squires transformation for vector fields
    for (symmTensor& r : R)
    {
        // (KSJ:Eq. 5)
        r.xx() = Foam::sqrt(r.xx());
        r.xy() /= r.xx();
        r.xz() /= r.xx();
        r.yy() = Foam::sqrt(r.yy() - sqr(r.xy()));
        r.yz() = (r.yz() - r.xy()*r.xz())/r.yy();
        r.zz() = Foam::sqrt(r.zz() - sqr(r.xz()) - sqr(r.yz()));
    }

    // Map transformed Reynolds stresses field onto patch for vector fields
    forAll(fld, i)
    {
        vector& u = fld[i];
        const symmTensor& r = R[i];

        // (KSJ:p. 658, item-e)
        u.z() = u.x()*r.xz() + u.y()*r.yz() + u.z()*r.zz();
        u.y() = u.x()*r.xy() + u.y()*r.yy();
        u.x() = u.x()*r.xx();
    }
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::mapMean
(
    scalarField& fld
) const
{
    const scalar t = this->db().time().timeOutputValue();

    fld += meanPtr_->value(t);
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::mapMean
(
    vectorField& fld
) const
{
    const scalar t = this->db().time().timeOutputValue();
    tmp<vectorField> tmean = meanPtr_->value(t);
    const vectorField& mean = tmean.cref();

    // Calculate flow-rate correction factor for vector fields (KCX:Eq. 8)
    const vector bulk
    (
        gSum(mean*this->patch().magSf())
       /(gSum(this->patch().magSf()) + ROOTVSMALL)
    );

    const scalar correct
    (
        gSum((bulk & patchNormal_)*this->patch().magSf())
       /gSum(mean & -this->patch().Sf())
    );

    // Map mean field onto patch for vector fields
    fld += mean;

    // Correct flow rate
    fld *= correct;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::turbulentDigitalFilterInletFvPatchField<Type>::
turbulentDigitalFilterInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    AMIPtr_(new faceAreaWeightAMI(true, false)),
    meanPtr_(nullptr),
    Rptr_(nullptr),
    curTimeIndex_(-1),
    patchNormal_(Zero),
    L_(p)
{}


template<class Type>
Foam::turbulentDigitalFilterInletFvPatchField<Type>::
turbulentDigitalFilterInletFvPatchField
(
    const turbulentDigitalFilterInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    AMIPtr_(ptf.AMIPtr_.clone()),
    meanPtr_(ptf.meanPtr_.clone(this->patch().patch())),
    Rptr_(ptf.Rptr_.clone(this->patch().patch())),
    curTimeIndex_(ptf.curTimeIndex_),
    patchNormal_(ptf.patchNormal_),
    L_(p, ptf.L_)
{}


template<class Type>
Foam::turbulentDigitalFilterInletFvPatchField<Type>::
turbulentDigitalFilterInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            true // flipNormals
        )
    ),
    meanPtr_
    (
        PatchFunction1<Type>::New
        (
            this->patch().patch(),
            "mean",
            dict
        )
    ),
    Rptr_
    (
        PatchFunction1<TypeR>::New
        (
            this->patch().patch(),
            "R",
            dict
        )
    ),
    curTimeIndex_(-1),
    patchNormal_(Zero),
    L_(p, dict)
{
    turbulence::IntegralScaleBox<Type>::debug = debug;

    // Check if varying or fixed time-step computation
    if (!L_.fsm() && this->db().time().isAdjustTimeStep())
    {
        WarningInFunction
            << "Varying time-step computations are not "
            << "supported by the digital filter method."
            << endl;
    }

    const scalar t = this->db().time().timeOutputValue();
    const Field<TypeR> R(Rptr_->value(t));

    turbulentDFSEMInletFvPatchVectorField::checkStresses(R);
}


template<class Type>
Foam::turbulentDigitalFilterInletFvPatchField<Type>::
turbulentDigitalFilterInletFvPatchField
(
    const turbulentDigitalFilterInletFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    AMIPtr_(ptf.AMIPtr_.clone()),
    meanPtr_(ptf.meanPtr_.clone(this->patch().patch())),
    Rptr_(ptf.Rptr_.clone(this->patch().patch())),
    curTimeIndex_(ptf.curTimeIndex_),
    patchNormal_(ptf.patchNormal_),
    L_(ptf.L_)
{}


template<class Type>
Foam::turbulentDigitalFilterInletFvPatchField<Type>::
turbulentDigitalFilterInletFvPatchField
(
    const turbulentDigitalFilterInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    AMIPtr_(ptf.AMIPtr_.clone()),
    meanPtr_(ptf.meanPtr_.clone(this->patch().patch())),
    Rptr_(ptf.Rptr_.clone(this->patch().patch())),
    curTimeIndex_(ptf.curTimeIndex_),
    patchNormal_(ptf.patchNormal_),
    L_(ptf.L_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);

    if (meanPtr_)
    {
        meanPtr_->autoMap(m);
    }
    if (Rptr_)
    {
        Rptr_->autoMap(m);
    }
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const auto& dfmptf =
        refCast<const turbulentDigitalFilterInletFvPatchField<Type>>(ptf);

    if (meanPtr_)
    {
        meanPtr_->rmap(dfmptf.meanPtr_(), addr);
    }
    if (Rptr_)
    {
        Rptr_->rmap(dfmptf.Rptr_(), addr);
    }
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& fld = *this;
        fld = Zero;

        mapL(fld);

        mapR(fld);

        mapMean(fld);

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::turbulentDigitalFilterInletFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);

    if (meanPtr_)
    {
        meanPtr_->writeData(os);
    }
    if (Rptr_)
    {
        Rptr_->writeData(os);
    }
    if (AMIPtr_)
    {
        AMIPtr_->write(os);
    }
    L_.write(os);

    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
