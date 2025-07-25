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

#include "uniformMixedFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    mixedFaPatchField<Type>(p, iF),
    refValueFunc_(nullptr),
    refGradFunc_(nullptr),
    valueFractionFunc_(nullptr)
{}


template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& fld
)
:
    mixedFaPatchField<Type>(p, iF, fld),
    refValueFunc_(nullptr),
    refGradFunc_(nullptr),
    valueFractionFunc_(nullptr)
{}


template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    // Bypass dict constructor, default initialise as zero-gradient
    mixedFaPatchField<Type>(p, iF, Foam::zero{}),
    refValueFunc_
    (
        Function1<Type>::NewIfPresent
        (
            /* p.patch(), */
            "uniformValue",
            dict,
            &iF.db()
        )
    ),
    refGradFunc_
    (
        Function1<Type>::NewIfPresent
        (
            // p.patch(),
            "uniformGradient",
            dict,
            &iF.db()
        )
    ),
    valueFractionFunc_(nullptr)
{
    faPatchFieldBase::readDict(dict);  // Consistent with a dict constructor

    if (refValueFunc_)
    {
        if (refGradFunc_)
        {
            // Both value + gradient: needs valueFraction
            valueFractionFunc_.reset
            (
                Function1<scalar>::New
                (
                    /* p.patch(), */
                    "uniformValueFraction",
                    dict,
                    &iF.db()
                )
            );
        }
    }
    else if (!refGradFunc_)
    {
        // Missing both value and gradient: FatalIOError
        FatalIOErrorInFunction(dict)
            << "For " << this->internalField().name() << " on "
            << this->patch().name() << nl
            << "Require either or both: uniformValue and uniformGradient"
            << " (possibly uniformValueFraction as well)" << nl
            << exit(FatalIOError);
    }

    // Use restart value if provided...
    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        // Evaluate to assign a value
        this->evaluate();
    }
}


template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const uniformMixedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    mixedFaPatchField<Type>(ptf, p, iF, mapper),
    refValueFunc_(ptf.refValueFunc_.clone(/*p.patch()*/)),
    refGradFunc_(ptf.refGradFunc_.clone(/*p.patch()*/)),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(/*p.patch()*/))
{}


template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const uniformMixedFaPatchField<Type>& ptf
)
:
    mixedFaPatchField<Type>(ptf),
    refValueFunc_(ptf.refValueFunc_.clone(/*this->patch().patch()*/)),
    refGradFunc_(ptf.refGradFunc_.clone(/*this->patch().patch()*/)),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(/*this->patch().patch()*/))
{}


template<class Type>
Foam::uniformMixedFaPatchField<Type>::uniformMixedFaPatchField
(
    const uniformMixedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    mixedFaPatchField<Type>(ptf, iF),
    refValueFunc_(ptf.refValueFunc_.clone(/*this->patch().patch()*/)),
    refGradFunc_(ptf.refGradFunc_.clone(/*this->patch().patch()*/)),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(/*this->patch().patch()*/))
{
    // Evaluate the profile if defined
    if (ptf.refValueFunc_ || ptf.refGradFunc_)
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformMixedFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();

    if (refValueFunc_)
    {
        this->refValue() = refValueFunc_->value(t);

        if (refGradFunc_)
        {
            // Both value + gradient: has valueFraction too
            this->valueFraction() = valueFractionFunc_->value(t);
        }
        else
        {
            // Has value only
            this->valueFraction() = 1;
        }
    }
    else
    {
        this->refValue() = Zero;
        this->valueFraction() = 0;
    }
    if (refGradFunc_)
    {
        this->refGrad() = refGradFunc_->value(t);
    }
    else
    {
        this->refGrad() = Zero;
    }

    // Missing both value and gradient is caught as an error in
    // dictionary constructor, but treated as zero-gradient here.

    mixedFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformMixedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);

    if (refValueFunc_)
    {
        refValueFunc_->writeData(os);
    }
    if (refGradFunc_)
    {
        refGradFunc_->writeData(os);
    }
    if (valueFractionFunc_)
    {
        valueFractionFunc_->writeData(os);
    }

    // For visualisation / restart
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
