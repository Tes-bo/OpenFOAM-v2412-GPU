/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "UniformValueField.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::PatchFunction1Types::UniformValueField<Type>::UniformValueField
(
    const polyPatch& pp,
    const word& redirectType,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    PatchFunction1<Type>(pp, entryName, dict, faceValues),
    uniformValuePtr_
    (
        Function1<Type>::New
        (
            entryName,
            dict,
            redirectType,
            patchFunction1Base::whichDb()  // mesh registry
        )
    )
{}


template<class Type>
Foam::PatchFunction1Types::UniformValueField<Type>::UniformValueField
(
    const UniformValueField<Type>& rhs,
    const polyPatch& pp
)
:
    PatchFunction1<Type>(rhs, pp),
    uniformValuePtr_(rhs.uniformValuePtr_.clone())
{
    if (uniformValuePtr_)
    {
        uniformValuePtr_->resetDb(patchFunction1Base::whichDb());
    }
}


template<class Type>
Foam::PatchFunction1Types::UniformValueField<Type>::UniformValueField
(
    const UniformValueField<Type>& rhs
)
:
    UniformValueField<Type>(rhs, rhs.patch())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::PatchFunction1Types::UniformValueField<Type>::autoMap
(
    const FieldMapper& mapper
)
{}


template<class Type>
void Foam::PatchFunction1Types::UniformValueField<Type>::rmap
(
    const PatchFunction1<Type>& pf1,
    const labelList& addr
)
{}


template<class Type>
void Foam::PatchFunction1Types::UniformValueField<Type>::writeData
(
    Ostream& os
) const
{
    PatchFunction1<Type>::writeData(os);

    uniformValuePtr_->writeData(os);
}


// ************************************************************************* //
