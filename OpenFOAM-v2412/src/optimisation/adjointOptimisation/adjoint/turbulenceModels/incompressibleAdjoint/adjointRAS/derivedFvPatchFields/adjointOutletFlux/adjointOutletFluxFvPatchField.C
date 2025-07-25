/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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

#include "adjointOutletFluxFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::adjointOutletFluxFvPatchField<Type>::adjointOutletFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
Foam::adjointOutletFluxFvPatchField<Type>::adjointOutletFluxFvPatchField
(
    const adjointOutletFluxFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
Foam::adjointOutletFluxFvPatchField<Type>::adjointOutletFluxFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, IOobjectOption::MUST_READ)
{}


template<class Type>
Foam::adjointOutletFluxFvPatchField<Type>::adjointOutletFluxFvPatchField
(
    const adjointOutletFluxFvPatchField<Type>& tppsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(tppsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::adjointOutletFluxFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    this->operator==(Field<Type>(this->patch().size(), Zero));

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::adjointOutletFluxFvPatchField<Type>::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::adjointOutletFluxFvPatchField<Type>::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::adjointOutletFluxFvPatchField<Type>::gradientBoundaryCoeffs() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
Foam::tmp<Foam::Field<Type>>
Foam::adjointOutletFluxFvPatchField<Type>::gradientInternalCoeffs() const
{
    return tmp<Field<Type>>::New(this->size(), Zero);
}


template<class Type>
void Foam::adjointOutletFluxFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
