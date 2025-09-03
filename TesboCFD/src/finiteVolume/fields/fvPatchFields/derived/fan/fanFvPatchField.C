/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "fanFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class Type>
const Foam::Enum<typename Foam::fanFvPatchField<Type>::operatingMode>
Foam::fanFvPatchField<Type>::operatingModeNames_
({
    { operatingMode::VELOCITY, "velocity" },
    { operatingMode::UNIFORM_VELOCITY, "uniformVelocity" },
    { operatingMode::VOL_FLOW_RATE, "volumeFlowRate" },
    { operatingMode::NON_DIMENSIONAL, "nonDimensional" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::fanFvPatchField<Type>::calcFanJump()
{
    if (this->cyclicPatch().owner())
    {
        this->jump_ = this->jumpTable_->value(this->db().time().value());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    uniformJumpFvPatchField<Type>(p, iF),
    operatingMode_(operatingMode::VELOCITY),
    phiName_("phi"),
    rhoName_("rho"),
    rpm_(nullptr),
    dm_(nullptr)
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    uniformJumpFvPatchField<Type>(p, iF, dict, false),  // needValue = false
    operatingMode_
    (
        operatingModeNames_.getOrDefault("mode", dict, operatingMode::VELOCITY)
    ),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho")),
    rpm_(nullptr),
    dm_(nullptr)
{
    // Backwards compatibility
    if (operatingMode_ == operatingMode::VELOCITY)
    {
        bool nonDimCompat = dict.getOrDefault("nonDimensional", false);
        if (nonDimCompat)
        {
            // Warn?
            operatingMode_ = operatingMode::NON_DIMENSIONAL;
        }

        bool uniformCompat = dict.getOrDefault("uniformJump", false);
        if (uniformCompat)
        {
            // Warn?
            operatingMode_ = operatingMode::UNIFORM_VELOCITY;
        }
    }


    // Note that we've not read jumpTable_ etc
    if (operatingMode_ == operatingMode::NON_DIMENSIONAL)
    {
        rpm_.reset(Function1<scalar>::New("rpm", dict, &this->db()));
        dm_.reset(Function1<scalar>::New("dm", dict, &this->db()));
    }

    if (this->cyclicPatch().owner())
    {
        this->jumpTable_ = Function1<Type>::New("jumpTable", dict, &this->db());
    }

    if (!this->readValueEntry(dict))
    {
        this->evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& rhs,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    uniformJumpFvPatchField<Type>(rhs, p, iF, mapper),
    operatingMode_(rhs.operatingMode_),
    phiName_(rhs.phiName_),
    rhoName_(rhs.rhoName_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& rhs
)
:
    uniformJumpFvPatchField<Type>(rhs),
    operatingMode_(rhs.operatingMode_),
    phiName_(rhs.phiName_),
    rhoName_(rhs.rhoName_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


template<class Type>
Foam::fanFvPatchField<Type>::fanFvPatchField
(
    const fanFvPatchField<Type>& rhs,
    const DimensionedField<Type, volMesh>& iF
)
:
    uniformJumpFvPatchField<Type>(rhs, iF),
    operatingMode_(rhs.operatingMode_),
    phiName_(rhs.phiName_),
    rhoName_(rhs.rhoName_),
    rpm_(rhs.rpm_.clone()),
    dm_(rhs.dm_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fanFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    calcFanJump();

    // Call fixedJump variant - uniformJump will overwrite the jump value
    fixedJumpFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::fanFvPatchField<Type>::write(Ostream& os) const
{
    uniformJumpFvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);

    os.writeEntry("mode", operatingModeNames_[operatingMode_]);

    if (operatingMode_ == operatingMode::NON_DIMENSIONAL)
    {
        rpm_->writeData(os);
        dm_->writeData(os);
    }
}


// ************************************************************************* //
