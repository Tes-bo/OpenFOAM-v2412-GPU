/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "uniformFixedValueFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF),
    refValueFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFaPatchField<Type>(p, iF, fld),
    refValueFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    refValueFunc_
    (
        Function1<Type>::New
        (
            /* p.patch(), */
            "uniformValue",
            dict
        )
    )
{
    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        // Evaluate to assign a value
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(p, iF),   // Don't map
    refValueFunc_(ptf.refValueFunc_.clone(/*p.patch()*/))
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf),
    refValueFunc_(ptf.refValueFunc_.clone(/*this->patch().patch()*/))
{}


template<class Type>
Foam::uniformFixedValueFaPatchField<Type>::uniformFixedValueFaPatchField
(
    const uniformFixedValueFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF),
    refValueFunc_(ptf.refValueFunc_.clone(/*this->patch().patch()*/))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    faPatchField<Type>::operator==(refValueFunc_->value(t));
    fixedValueFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedValueFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    if (refValueFunc_)
    {
        refValueFunc_->writeData(os);
    }
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
