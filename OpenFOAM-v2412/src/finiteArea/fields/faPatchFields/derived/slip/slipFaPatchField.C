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

#include "slipFaPatchField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::slipFaPatchField<Type>::slipFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    basicSymmetryFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::slipFaPatchField<Type>::slipFaPatchField
(
    const slipFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    basicSymmetryFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::slipFaPatchField<Type>::slipFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    basicSymmetryFaPatchField<Type>(p, iF, dict)
{}


template<class Type>
Foam::slipFaPatchField<Type>::slipFaPatchField
(
    const slipFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    basicSymmetryFaPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::slipFaPatchField<Type>::slipFaPatchField
(
    const slipFaPatchField<Type>& ptf
)
:
    basicSymmetryFaPatchField<Type>(ptf)
{}


// ************************************************************************* //
