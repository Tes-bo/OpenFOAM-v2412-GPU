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

#include "fixedValueOutflowFaPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict, IOobjectOption::MUST_READ)
{}


template<class Type>
Foam::fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


template<class Type>
Foam::fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::fixedValueOutflowFaPatchField<Type>::fixedValueOutflowFaPatchField
(
    const fixedValueOutflowFaPatchField<Type>& ptf
)
:
    faPatchField<Type>(ptf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueOutflowFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>& weights
) const
{
    return pTraits<Type>::one*weights;
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueOutflowFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>& weights
) const
{
    return (1.0 - weights)*(*this);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueOutflowFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -Type(pTraits<Type>::one)*this->patch().deltaCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::fixedValueOutflowFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void Foam::fixedValueOutflowFaPatchField<Type>::write(Ostream& os) const
{
    faPatchField<Type>::write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
