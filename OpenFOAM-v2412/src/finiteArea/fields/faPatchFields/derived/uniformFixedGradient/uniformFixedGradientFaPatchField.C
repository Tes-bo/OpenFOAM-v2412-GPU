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

#include "uniformFixedGradientFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedGradientFaPatchField<Type>(p, iF),
    refGradFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& fld
)
:
    fixedGradientFaPatchField<Type>(p, iF, fld),
    refGradFunc_(nullptr)
{}


template<class Type>
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFaPatchField<Type>(p, iF),  // Bypass dictionary constructor
    refGradFunc_
    (
        Function1<Type>::New
        (
            /* p.patch(), */
            "uniformGradient",
            dict,
            &iF.db()
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
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const uniformFixedGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<Type>(ptf, p, iF, mapper),
    refGradFunc_(ptf.refGradFunc_.clone(/*p.patch()*/))
{}


template<class Type>
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const uniformFixedGradientFaPatchField<Type>& ptf
)
:
    fixedGradientFaPatchField<Type>(ptf),
    refGradFunc_(ptf.refGradFunc_.clone(/*p.patch()*/))
{}


template<class Type>
Foam::uniformFixedGradientFaPatchField<Type>::uniformFixedGradientFaPatchField
(
    const uniformFixedGradientFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedGradientFaPatchField<Type>(ptf, iF),
    refGradFunc_(ptf.refGradFunc_.clone(/*this->patch().patch()*/))
{
    // Evaluate the profile if defined
    if (refGradFunc_)
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformFixedGradientFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    // Extra safety on the evaluation
    if (refGradFunc_)
    {
        this->gradient() = refGradFunc_->value(t);
    }
    else
    {
        this->gradient() = Zero;
    }

    fixedGradientFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformFixedGradientFaPatchField<Type>::write(Ostream& os) const
{
    fixedGradientFaPatchField<Type>::write(os);
    if (refGradFunc_)
    {
        refGradFunc_->writeData(os);
    }
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
