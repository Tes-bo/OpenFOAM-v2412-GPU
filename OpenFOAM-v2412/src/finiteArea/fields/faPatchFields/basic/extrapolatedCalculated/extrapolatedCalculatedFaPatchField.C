/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "extrapolatedCalculatedFaPatchField.H"
#include "faPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::extrapolatedCalculatedFaPatchField<Type>::
extrapolatedCalculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    calculatedFaPatchField<Type>(p, iF)
{}


template<class Type>
Foam::extrapolatedCalculatedFaPatchField<Type>::
extrapolatedCalculatedFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    calculatedFaPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ)
{
    faPatchField<Type>::extrapolateInternal();  // Zero-gradient patch values
}


template<class Type>
Foam::extrapolatedCalculatedFaPatchField<Type>::
extrapolatedCalculatedFaPatchField
(
    const extrapolatedCalculatedFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    calculatedFaPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::extrapolatedCalculatedFaPatchField<Type>::
extrapolatedCalculatedFaPatchField
(
    const extrapolatedCalculatedFaPatchField<Type>& ptf
)
:
    calculatedFaPatchField<Type>(ptf)
{}


template<class Type>
Foam::extrapolatedCalculatedFaPatchField<Type>::
extrapolatedCalculatedFaPatchField
(
    const extrapolatedCalculatedFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    calculatedFaPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::extrapolatedCalculatedFaPatchField<Type>::evaluate
(
    const Pstream::commsTypes
)
{
    if (!this->updated())
    {
        this->updateCoeffs();
    }

    faPatchField<Type>::extrapolateInternal();  // Zero-gradient patch values
    calculatedFaPatchField<Type>::evaluate();
}


// ************************************************************************* //
