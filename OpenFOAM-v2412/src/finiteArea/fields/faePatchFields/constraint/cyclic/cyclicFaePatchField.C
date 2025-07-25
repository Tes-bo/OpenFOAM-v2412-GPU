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

#include "cyclicFaePatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::cyclicFaePatchField<Type>::cyclicFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    coupledFaePatchField<Type>(p, iF),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{}


template<class Type>
Foam::cyclicFaePatchField<Type>::cyclicFaePatchField
(
    const cyclicFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    coupledFaePatchField<Type>(ptf, p, iF, mapper),
    cyclicPatch_(refCast<const cyclicFaPatch>(p))
{
    if (!isType<cyclicFaPatch>(this->patch()))
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
Foam::cyclicFaePatchField<Type>::cyclicFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    coupledFaePatchField<Type>(p, iF, dict),
    cyclicPatch_(refCast<const cyclicFaPatch>(p, dict))
{
    if (!isType<cyclicFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not cyclic type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::cyclicFaePatchField<Type>::cyclicFaePatchField
(
    const cyclicFaePatchField<Type>& ptf
)
:
    coupledFaePatchField<Type>(ptf),
    cyclicPatch_(ptf.cyclicPatch_)
{}


template<class Type>
Foam::cyclicFaePatchField<Type>::cyclicFaePatchField
(
    const cyclicFaePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    coupledFaePatchField<Type>(ptf, iF),
    cyclicPatch_(ptf.cyclicPatch_)
{}


// ************************************************************************* //
