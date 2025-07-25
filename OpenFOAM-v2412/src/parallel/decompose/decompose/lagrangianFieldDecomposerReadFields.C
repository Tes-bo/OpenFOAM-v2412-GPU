/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "lagrangianFieldDecomposer.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::lagrangianFieldDecomposer::readFields
(
    const label cloudi,
    const IOobjectList& lagrangianObjects,
    PtrList<PtrList<IOField<Type>>>& lagrangianFields
)
{
    // List of lagrangian field objects
    const UPtrList<const IOobject> fieldObjects
    (
        lagrangianObjects.csorted<IOField<Type>>()
    );

    auto& cloudFields =
        lagrangianFields.emplace_set(cloudi, fieldObjects.size());

    forAll(fieldObjects, fieldi)
    {
        cloudFields.emplace_set(fieldi, fieldObjects[fieldi]);
    }
}


template<class Type>
void Foam::lagrangianFieldDecomposer::readFieldFields
(
    const label cloudi,
    const IOobjectList& lagrangianObjects,
    PtrList<PtrList<CompactIOField<Field<Type>, Type>>>& lagrangianFields
)
{
    // List of lagrangian field objects
    UPtrList<const IOobject> fieldObjects
    (
        lagrangianObjects.cobjects<IOField<Field<Type>>>()
    );

    fieldObjects.push_back
    (
        lagrangianObjects.cobjects<CompactIOField<Field<Type>, Type>>()
    );

    Foam::sort(fieldObjects, nameOp<IOobject>());


    auto& cloudFields =
        lagrangianFields.emplace_set(cloudi, fieldObjects.size());

    forAll(fieldObjects, fieldi)
    {
        cloudFields.emplace_set(fieldi, fieldObjects[fieldi]);
    }
}


// ************************************************************************* //
