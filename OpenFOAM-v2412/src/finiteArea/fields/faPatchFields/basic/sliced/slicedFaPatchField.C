/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "slicedFaPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const Field<Type>& completeOrBoundaryField,
    const bool isBoundaryOnly
)
:
    faPatchField<Type>(p, iF, Field<Type>())
{
    if (isBoundaryOnly)
    {
        // Set to a slice of the boundary field
        UList<Type>::shallowCopy(p.boundarySlice(completeOrBoundaryField));
    }
    else
    {
        // Set to a slice of the complete field
        UList<Type>::shallowCopy(p.patchSlice(completeOrBoundaryField));
    }
}


template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF, Field<Type>())
{}


template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF)  // bypass dictionary constructor
{
    faPatchFieldBase::readDict(dict);
    // Read "value" if present...

    NotImplemented;
}


template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const slicedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{
    NotImplemented;
}


template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const slicedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf.patch(), iF, Field<Type>())
{
    // Transfer the slice from the argument
    UList<Type>::shallowCopy(ptf);
}


template<class Type>
Foam::slicedFaPatchField<Type>::slicedFaPatchField
(
    const slicedFaPatchField<Type>& ptf
)
:
    slicedFaPatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::slicedFaPatchField<Type>::~slicedFaPatchField()
{
    // Set to nullptr to avoid deletion of underlying field
    UList<Type>::shallowCopy(nullptr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::slicedFaPatchField<Type>::snGrad() const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::patchInternalField() const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::slicedFaPatchField<Type>::patchInternalField(Field<Type>&) const
{
    NotImplemented;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::patchNeighbourField
(
    const Field<Type>& iField
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::patchNeighbourField() const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::gradientInternalCoeffs() const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::slicedFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    NotImplemented;
    return nullptr;
}


template<class Type>
void Foam::slicedFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
