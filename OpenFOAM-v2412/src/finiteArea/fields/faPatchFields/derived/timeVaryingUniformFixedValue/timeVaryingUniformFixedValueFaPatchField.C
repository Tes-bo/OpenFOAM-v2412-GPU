/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "timeVaryingUniformFixedValueFaPatchField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(p, iF),
    timeSeries_()
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFaPatchField<Type>(p, iF, dict, IOobjectOption::NO_READ),
    timeSeries_(dict)
{
    if (!this->readValueEntry(dict))
    {
        // Ensure field has reasonable initial values
        this->extrapolateInternal();

        updateCoeffs();
    }

    DeprecatedInFunction(2212)
        << "Use uniformFixedValue instead." << nl
        << "    This boundary condition will be removed in the future" << endl;
}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf,
    const faPatch& p,
    const DimensionedField<Type, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedValueFaPatchField<Type>(ptf, p, iF, mapper),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf
)
:
    fixedValueFaPatchField<Type>(ptf),
    timeSeries_(ptf.timeSeries_)
{}


template<class Type>
Foam::timeVaryingUniformFixedValueFaPatchField<Type>::
timeVaryingUniformFixedValueFaPatchField
(
    const timeVaryingUniformFixedValueFaPatchField<Type>& ptf,
    const DimensionedField<Type, areaMesh>& iF
)
:
    fixedValueFaPatchField<Type>(ptf, iF),
    timeSeries_(ptf.timeSeries_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::timeVaryingUniformFixedValueFaPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    faPatchField<Type>::operator==
    (
        timeSeries_(this->db().time().timeOutputValue())
    );
    fixedValueFaPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::timeVaryingUniformFixedValueFaPatchField<Type>::write
(
    Ostream& os
) const
{
    faPatchField<Type>::write(os);
    timeSeries_.write(os);
    faPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
