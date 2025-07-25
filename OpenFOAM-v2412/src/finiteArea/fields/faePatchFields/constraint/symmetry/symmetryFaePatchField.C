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

#include "symmetryFaePatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::symmetryFaePatchField<Type>::symmetryFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(p, iF)
{}


template<class Type>
Foam::symmetryFaePatchField<Type>::symmetryFaePatchField
(
    const symmetryFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchField<Type>(ptf, p, iF, mapper)
{
    if (!isType<symmetryFaPatch>(this->patch()))
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
Foam::symmetryFaePatchField<Type>::symmetryFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict
)
:
    faePatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{
    if (!isType<symmetryFaPatch>(p))
    {
        FatalIOErrorInFunction(dict)
            << "patch " << this->patch().index() << " not symmetry type. "
            << "Patch type = " << p.type()
            << exit(FatalIOError);
    }
}


template<class Type>
Foam::symmetryFaePatchField<Type>::symmetryFaePatchField
(
    const symmetryFaePatchField<Type>& ptf
)
:
    faePatchField<Type>(ptf)
{}


template<class Type>
Foam::symmetryFaePatchField<Type>::symmetryFaePatchField
(
    const symmetryFaePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(ptf, iF)
{}


// ************************************************************************* //
