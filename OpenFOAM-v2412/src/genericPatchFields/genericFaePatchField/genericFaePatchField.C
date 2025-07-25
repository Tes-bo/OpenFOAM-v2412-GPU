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

#include "genericFaePatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::genericFaePatchField<Type>::genericFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    parent_bctype(p, iF)
{
    FatalErrorInFunction
        << "Trying to construct generic patchField on patch "
        << this->patch().name()
        << " of field " << this->internalField().name()
        << abort(FatalError);
}


template<class Type>
Foam::genericFaePatchField<Type>::genericFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict),
    genericPatchFieldBase(dict)
{
    const label patchSize = this->size();
    const word& patchName = this->patch().name();
    const IOobject& io = this->internalField();

    if (!dict.findEntry("value", keyType::LITERAL))
    {
        reportMissingEntry("value", patchName, io);
    }

    // Handle "value" separately
    processGeneric(patchSize, patchName, io, true);
}


template<class Type>
Foam::genericFaePatchField<Type>::genericFaePatchField
(
    const genericFaePatchField<Type>& rhs,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper),
    genericPatchFieldBase(zero{}, rhs)
{
    this->mapGeneric(rhs, mapper);
}


template<class Type>
Foam::genericFaePatchField<Type>::genericFaePatchField
(
    const genericFaePatchField<Type>& rhs,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    parent_bctype(rhs, iF),
    genericPatchFieldBase(rhs)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::genericFaePatchField<Type>::write(Ostream& os) const
{
    // Handle "value" separately
    genericPatchFieldBase::writeGeneric(os, true);
    faePatchField<Type>::writeValueEntry(os);
}


template<class Type>
void Foam::genericFaePatchField<Type>::autoMap
(
    const faPatchFieldMapper& m
)
{
    parent_bctype::autoMap(m);
    this->autoMapGeneric(m);
}


template<class Type>
void Foam::genericFaePatchField<Type>::rmap
(
    const faePatchField<Type>& rhs,
    const labelList& addr
)
{
    parent_bctype::rmap(rhs, addr);

    const auto* base = isA<genericPatchFieldBase>(rhs);
    if (base)
    {
        this->rmapGeneric(*base, addr);
    }
}


// ************************************************************************* //
