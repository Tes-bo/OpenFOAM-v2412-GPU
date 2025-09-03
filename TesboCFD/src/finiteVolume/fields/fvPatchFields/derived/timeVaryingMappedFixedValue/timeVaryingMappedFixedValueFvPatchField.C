/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2018 OpenCFD Ltd.
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

#include "timeVaryingMappedFixedValueFvPatchField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    uniformValue_()
{}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            p.patch(),
            "uniformValue",
            dict,
            iF.name(),          // field table name
            true                // face values
        )
    )
{
    if (!this->readValueEntry(dict))
    {
        // Note: we use evaluate() here to trigger updateCoeffs followed
        //       by re-setting of fvatchfield::updated_ flag. This is
        //       so if first use is in the next time step it retriggers
        //       a new update.
        this->evaluate(Pstream::commsTypes::buffered);
    }
}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            p.patch()
        )
    )
{}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            this->patch().patch()
        )
    )
{}


template<class Type>
Foam::timeVaryingMappedFixedValueFvPatchField<Type>::
timeVaryingMappedFixedValueFvPatchField
(
    const timeVaryingMappedFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    uniformValue_
    (
        new PatchFunction1Types::MappedFile<Type>
        (
            ptf.uniformValue_(),
            this->patch().patch()
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    uniformValue_().autoMap(m);
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const timeVaryingMappedFixedValueFvPatchField<Type>& tiptf =
        refCast<const timeVaryingMappedFixedValueFvPatchField<Type>>(ptf);

    uniformValue_().rmap(tiptf.uniformValue_(), addr);
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }


    const scalar t = this->db().time().timeOutputValue();
    fvPatchField<Type>::operator==(uniformValue_->value(t));

    if (debug)
    {
        Pout<< "updateCoeffs : set fixedValue to min:" << gMin(*this)
            << " max:" << gMax(*this)
            << " avg:" << gAverage(*this) << endl;
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingMappedFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    uniformValue_->writeData(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
