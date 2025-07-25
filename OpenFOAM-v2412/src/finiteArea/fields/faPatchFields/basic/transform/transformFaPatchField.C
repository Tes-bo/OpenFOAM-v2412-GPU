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

#include "transformFaPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::transformFaPatchField<Type>::transformFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
   faPatchField<Type>(p, iF)
{}


template<class Type>
Foam::transformFaPatchField<Type>::transformFaPatchField
(
    const transformFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    faPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::transformFaPatchField<Type>::transformFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    faPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{}


template<class Type>
Foam::transformFaPatchField<Type>::transformFaPatchField
(
    const transformFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    faPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFaPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return pTraits<Type>::one - snGradTransformDiag();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFaPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return
        *this
      - cmptMultiply
        (
            valueInternalCoeffs(this->patch().weights()),
            this->patchInternalField()
        );
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFaPatchField<Type>::gradientInternalCoeffs() const
{
    return -this->patch().deltaCoeffs()*snGradTransformDiag();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::transformFaPatchField<Type>::gradientBoundaryCoeffs() const
{
    return
        snGrad()
      - cmptMultiply(gradientInternalCoeffs(), this->patchInternalField());
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::transformFaPatchField<Type>::operator=
(
    const faPatchField<Type>& ptf
)
{
    this->evaluate();
}


// ************************************************************************* //
