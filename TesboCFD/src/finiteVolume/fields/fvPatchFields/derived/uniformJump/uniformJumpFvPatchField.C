/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

#include "uniformJumpFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(p, iF),
    jumpTable_()
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedJumpFvPatchField<Type>(ptf, p, iF, mapper),
    jumpTable_(ptf.jumpTable_.clone())
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict,
    const bool needValue
)
:
    fixedJumpFvPatchField<Type>(p, iF, dict, false),  // needValue = false
    jumpTable_()
{
    if (needValue)
    {
        if (this->cyclicPatch().owner())
        {
            jumpTable_ = Function1<Type>::New("jumpTable", dict, &this->db());
        }

        if (!this->readValueEntry(dict))
        {
            this->evaluate(Pstream::commsTypes::buffered);
        }
    }
}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf
)
:
    fixedJumpFvPatchField<Type>(ptf),
    jumpTable_(ptf.jumpTable_.clone())
{}


template<class Type>
Foam::uniformJumpFvPatchField<Type>::uniformJumpFvPatchField
(
    const uniformJumpFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedJumpFvPatchField<Type>(ptf, iF),
    jumpTable_(ptf.jumpTable_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformJumpFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (this->cyclicPatch().owner())
    {
        this->setJump(jumpTable_->value(this->db().time().value()));
    }

    fixedJumpFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformJumpFvPatchField<Type>::write(Ostream& os) const
{
    fixedJumpFvPatchField<Type>::write(os);

    if (this->cyclicPatch().owner())
    {
        jumpTable_->writeData(os);
    }
}


// ************************************************************************* //
