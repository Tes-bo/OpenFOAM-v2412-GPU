/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "zeroFixedValuePointPatchField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * * //

template<class Type>
Foam::zeroFixedValuePointPatchField<Type>::zeroFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    // Field is zero
    fixedValuePointPatchField<Type>(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroFixedValuePointPatchField<Type>::zeroFixedValuePointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    // Field is zero
    fixedValuePointPatchField<Type>(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroFixedValuePointPatchField<Type>::zeroFixedValuePointPatchField
(
    const zeroFixedValuePointPatchField<Type>& ptf,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    // Field is zero : no mapping needed
    fixedValuePointPatchField<Type>(p, iF, Type(Zero))
{}


template<class Type>
Foam::zeroFixedValuePointPatchField<Type>::zeroFixedValuePointPatchField
(
    const zeroFixedValuePointPatchField<Type>& ptf
)
:
    fixedValuePointPatchField<Type>(ptf)
{
    // Field is zero. For safety, re-assign
    valuePointPatchField<Type>::operator=(Type(Zero));
}


template<class Type>
Foam::zeroFixedValuePointPatchField<Type>::zeroFixedValuePointPatchField
(
    const zeroFixedValuePointPatchField<Type>& ptf,
    const DimensionedField<Type, pointMesh>& iF
)
:
    fixedValuePointPatchField<Type>(ptf, iF)
{
    // Field is zero. For safety, re-assign
    valuePointPatchField<Type>::operator=(Type(Zero));
}


// ************************************************************************* //
