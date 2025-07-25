/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "emptyFaePatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::emptyFaePatchField<Type>::emptyFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{}


template<class Type>
Foam::emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>&,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper&
)
:
    faePatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    if (!isType<emptyFaPatch>(this->patch()))
    {
        FatalErrorInFunction
            << "Field type does not correspond to patch type for patch "
            << this->patch().index() << "." << endl
            << "Field type: " << typeName << endl
            << "Patch type: " << this->patch().type()
            << exit(FatalError);
    }
}


template<class Type>
Foam::emptyFaePatchField<Type>::emptyFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, Field<Type>())  // zero-sized patch field
{
    // Empty means empty, so no patchType override
    // with faePatchFieldBase::readDict(dict);

    if (!isType<emptyFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not empty type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(ptf.patch(), iF, Field<Type>())  // zero-sized
{}


template<class Type>
Foam::emptyFaePatchField<Type>::emptyFaePatchField
(
    const emptyFaePatchField<Type>& ptf
)
:
    emptyFaePatchField<Type>(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::emptyFaePatchField<Type>::write(Ostream& os) const
{
    faePatchField<Type>::write(os);
    // Never write "value"
}


// ************************************************************************* //
