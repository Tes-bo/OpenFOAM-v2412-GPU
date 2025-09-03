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

#include "uniformMixedFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    refValueFunc_(nullptr),
    refGradFunc_(nullptr),
    valueFractionFunc_(nullptr)
{}


template<class Type>
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    mixedFvPatchField<Type>(p, iF, fld),
    refValueFunc_(nullptr),
    refGradFunc_(nullptr),
    valueFractionFunc_(nullptr)
{}


template<class Type>
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    // Bypass dict constructor, default initialise as zero-gradient
    mixedFvPatchField<Type>(p, iF, Foam::zero{}),
    refValueFunc_
    (
        PatchFunction1<Type>::NewIfPresent(p.patch(), "uniformValue", dict)
    ),
    refGradFunc_
    (
        PatchFunction1<Type>::NewIfPresent(p.patch(), "uniformGradient", dict)
    ),
    valueFractionFunc_(nullptr)
{
    fvPatchFieldBase::readDict(dict);  // Consistent with a dict constructor

    if (refValueFunc_)
    {
        if (refGradFunc_)
        {
            // Both value + gradient: needs valueFraction
            valueFractionFunc_.reset
            (
                PatchFunction1<scalar>::New
                (
                    p.patch(),
                    "uniformValueFraction",
                    dict
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
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const uniformMixedFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    refValueFunc_(ptf.refValueFunc_.clone(p.patch())),
    refGradFunc_(ptf.refGradFunc_.clone(p.patch())),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(p.patch()))
{}


template<class Type>
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const uniformMixedFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    refValueFunc_(ptf.refValueFunc_.clone(this->patch().patch())),
    refGradFunc_(ptf.refGradFunc_.clone(this->patch().patch())),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(this->patch().patch()))
{}


template<class Type>
Foam::uniformMixedFvPatchField<Type>::uniformMixedFvPatchField
(
    const uniformMixedFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    refValueFunc_(ptf.refValueFunc_.clone(this->patch().patch())),
    refGradFunc_(ptf.refGradFunc_.clone(this->patch().patch())),
    valueFractionFunc_(ptf.valueFractionFunc_.clone(this->patch().patch()))
{
    // Evaluate the profile if defined
    if (ptf.refValueFunc_ || ptf.refGradFunc_)
    {
        this->evaluate();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::uniformMixedFvPatchField<Type>::updateCoeffs()
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

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::uniformMixedFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

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
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
