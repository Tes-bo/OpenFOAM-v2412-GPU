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

#include "coupledFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(p, iF)
{}


template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const Field<Type>& f
)
:
    faePatchField<Type>(p, iF, f)
{}


template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faePatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const faPatch& p,
    const DimensionedField<Type, edgeMesh>& iF,
    const dictionary& dict,
    IOobjectOption::readOption requireValue
)
:
    faePatchField<Type>(p, iF, dict, requireValue)
{}


template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf
)
:
    faePatchField<Type>(ptf)
{}


template<class Type>
Foam::coupledFaePatchField<Type>::coupledFaePatchField
(
    const coupledFaePatchField<Type>& ptf,
    const DimensionedField<Type, edgeMesh>& iF
)
:
    faePatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::coupledFaePatchField<Type>::write(Ostream& os) const
{
    faePatchField<Type>::write(os);
    faePatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
