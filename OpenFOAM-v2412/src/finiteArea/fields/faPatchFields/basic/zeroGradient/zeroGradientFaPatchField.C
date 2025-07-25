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

#include "zeroGradientFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const zeroGradientFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{
    faPatchField<Type>::extrapolateInternal();  // Zero-gradient patch values
}


template<class Type>
Foam::zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const zeroGradientFaPatchField& zgpf
)
:
    faPatchField<Type>(zgpf)
{}


template<class Type>
Foam::zeroGradientFaPatchField<Type>::zeroGradientFaPatchField
(
    const zeroGradientFaPatchField& zgpf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(zgpf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::zeroGradientFaPatchField<Type>::evaluate(const Pstream::commsTypes)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    faPatchField<Type>::extrapolateInternal();  // Zero-gradient patch values
    faPatchField<Type>::evaluate();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroGradientFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>::New(this->size(), pTraits<Type>::one);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroGradientFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroGradientFaPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::zeroGradientFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


// ************************************************************************* //
