/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "outletInletFvPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::outletInletFvPatchField<Type>::outletInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_("phi")
{
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0;
}


template<class Type>
Foam::outletInletFvPatchField<Type>::outletInletFvPatchField
(
    const outletInletFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::outletInletFvPatchField<Type>::outletInletFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<Type>(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    fvPatchFieldBase::readDict(dict);

    // Require outletValue (MUST_READ)
    this->refValue().assign("outletValue", dict, p.size());
    this->refGrad() = Zero;
    this->valueFraction() = 0;

    if (!this->readValueEntry(dict))
    {
        fvPatchField<Type>::extrapolateInternal();
    }
}


template<class Type>
Foam::outletInletFvPatchField<Type>::outletInletFvPatchField
(
    const outletInletFvPatchField<Type>& ptf
)
:
    mixedFvPatchField<Type>(ptf),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::outletInletFvPatchField<Type>::outletInletFvPatchField
(
    const outletInletFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    mixedFvPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::outletInletFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<surfaceScalarField>(phiName_);

    this->valueFraction() = pos0(phip);

    mixedFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::outletInletFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    this->refValue().writeEntry("outletValue", os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
