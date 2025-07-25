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

#include "outletInletFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::outletInletFaPatchField<Type>::outletInletFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    mixedFaPatchField<Type>(p, iF),
    phiName_("phi")
{
    this->refValue() = *this;
    this->refGrad() = Zero;
    this->valueFraction() = 0;
}


template<class Type>
Foam::outletInletFaPatchField<Type>::outletInletFaPatchField
(
    const outletInletFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    mixedFaPatchField<Type>(ptf, p, iF, mapper),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::outletInletFaPatchField<Type>::outletInletFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    mixedFaPatchField<Type>(p, iF),
    phiName_(dict.getOrDefault<word>("phi", "phi"))
{
    faPatchFieldBase::readDict(dict);

    // Require outletValue (MUST_READ)
    this->refValue().assign("outletValue", dict, p.size());
    this->refGrad() = Zero;
    this->valueFraction() = 0;

    if (!this->readValueEntry(dict))
    {
        faPatchField<Type>::extrapolateInternal();
    }
}


template<class Type>
Foam::outletInletFaPatchField<Type>::outletInletFaPatchField
(
    const outletInletFaPatchField<Type>& ptf
)
:
    mixedFaPatchField<Type>(ptf),
    phiName_(ptf.phiName_)
{}


template<class Type>
Foam::outletInletFaPatchField<Type>::outletInletFaPatchField
(
    const outletInletFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    mixedFaPatchField<Type>(ptf, iF),
    phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::outletInletFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const Field<scalar>& phip =
        this->patch().template lookupPatchField<edgeScalarField>(phiName_);

    this->valueFraction() = pos0(phip);

    mixedFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::outletInletFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    this->refValue().writeEntry("outletValue", os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
