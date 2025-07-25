/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "genericPointPatchField.H"
#include "pointPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF
)
:
    parent_bctype(p, iF)
{
    FatalErrorInFunction
        << "Trying to construct generic patchField on patch "
        << this->patch().name()
        << " of field " << this->internalField().name() << nl
        << abort(FatalError);
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    genericPatchFieldBase(dict)
{
    const label patchSize = this->size();
    const word& patchName = this->patch().name();
    const IOobject& io = this->internalField();

    // No separate "value"
    processGeneric(patchSize, patchName, io, false);
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& rhs,
    const pointPatch& p,
    const DimensionedField<Type, pointMesh>& iF,
    const pointPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    genericPatchFieldBase(zero{}, rhs)
{
    this->mapGeneric(rhs, mapper);
}


template<class Type>
Foam::genericPointPatchField<Type>::genericPointPatchField
(
    const genericPointPatchField<Type>& rhs,
    const DimensionedField<Type, pointMesh>& iF
)
:
    parent_bctype(rhs, iF),
    genericPatchFieldBase(rhs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericPointPatchField<Type>::write(Ostream& os) const
{
    // No separate treatment for "value"
    genericPatchFieldBase::writeGeneric(os, false);
}


template<class Type>
void Foam::genericPointPatchField<Type>::autoMap
(
    const pointPatchFieldMapper& m
)
{
    this->autoMapGeneric(m);
}


template<class Type>
void Foam::genericPointPatchField<Type>::rmap
(
    const pointPatchField<Type>& rhs,
    const labelList& addr
)
{
    const auto* base = isA<genericPatchFieldBase>(rhs);
    if (base)
    {
        this->rmapGeneric(*base, addr);
    }
}


// ************************************************************************* //
