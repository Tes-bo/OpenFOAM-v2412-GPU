/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "acousticWaveTransmissiveFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(p, iF),
    advectiveU_(0)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    advectiveFvPatchField<Type>(ptf, p, iF, mapper),
    advectiveU_(ptf.advectiveU_)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    advectiveFvPatchField<Type>(p, iF, dict),
    advectiveU_(dict.get<scalar>("advectiveSpeed"))
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptpsf
)
:
    advectiveFvPatchField<Type>(ptpsf),
    advectiveU_(ptpsf.advectiveU_)
{}


template<class Type>
Foam::acousticWaveTransmissiveFvPatchField<Type>::acousticWaveTransmissiveFvPatchField
(
    const acousticWaveTransmissiveFvPatchField& ptpsf,
    const DimensionedField<Type, volMesh>& iF
)
:
    advectiveFvPatchField<Type>(ptpsf, iF),
    advectiveU_(ptpsf.advectiveU_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::scalarField>
Foam::acousticWaveTransmissiveFvPatchField<Type>::advectionSpeed() const
{
    return tmp<scalarField>::New(this->size(), advectiveU_);
}


template<class Type>
void Foam::acousticWaveTransmissiveFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeEntry("advectiveSpeed", advectiveU_);
    fvPatchField<Type>::writeValueEntry(os);
}


// ************************************************************************* //
