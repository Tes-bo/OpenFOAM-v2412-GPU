/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "surfMesh.H"
#include "surfFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class GeoMeshType>
Foam::DimensionedField<Type, GeoMeshType>&
Foam::surfMesh::newField
(
    const word& fieldName,
    const dimensionSet& dims
)
{
    typedef DimensionedField<Type, GeoMeshType> fieldType;

    const objectRegistry& fieldDb = *this;

    auto* fldptr = fieldDb.getObjectPtr<fieldType>(fieldName);

    if (fldptr)
    {
        fldptr->dimensions().reset(dims);  // Dimensions may have changed
        fldptr->field() = Foam::zero{};
    }
    else
    {
        fldptr = new fieldType
        (
            fieldDb.newIOobject
            (
                fieldName,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            *this,
            Foam::zero{},
            dims
        );

        regIOobject::store(fldptr);
    }

    return *fldptr;
}


template<class Type, class GeoMeshType>
Foam::DimensionedField<Type, GeoMeshType>&
Foam::surfMesh::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    const Field<Type>& values
)
{
    typedef DimensionedField<Type, GeoMeshType> fieldType;

    const objectRegistry& fieldDb = *this;

    auto* fldptr = fieldDb.getObjectPtr<fieldType>(fieldName);

    if (fldptr)
    {
        fldptr->dimensions().reset(dims);  // Dimensions may have changed
        fldptr->field() = values;
    }
    else
    {
        fldptr = new fieldType
        (
            fieldDb.newIOobject
            (
                fieldName,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            *this,
            dims,
            values
        );

        regIOobject::store(fldptr);
    }

    return *fldptr;
}


template<class Type, class GeoMeshType>
Foam::DimensionedField<Type, GeoMeshType>&
Foam::surfMesh::storeField
(
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values
)
{
    typedef DimensionedField<Type, GeoMeshType> fieldType;

    const objectRegistry& fieldDb = *this;

    auto* fldptr = fieldDb.getObjectPtr<fieldType>(fieldName);

    if (fldptr)
    {
        fldptr->dimensions().reset(dims);  // Dimensions may have changed
        fldptr->field() = std::move(values);
    }
    else
    {
        fldptr = new fieldType
        (
            fieldDb.newIOobject
            (
                fieldName,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            ),
            *this,
            dims,
            std::move(values)
        );

        regIOobject::store(fldptr);
    }

    return *fldptr;
}


// ************************************************************************* //
