/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "IOField.H"
#include "CompactIOField.H"
#include "Time.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::IOField<Type>>
Foam::lagrangianReconstructor::reconstructField
(
    const word& cloudName,
    const word& fieldName
)
{
    // Construct empty field on mesh
    auto tfield = tmp<IOField<Type>>::New
    (
        IOobject
        (
            fieldName,
            mesh_.time().timeName(),
            cloud::prefix/cloudName,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Field<Type>()
    );
    auto& field = tfield.ref();

    for (const fvMesh& localMesh : procMeshes_)
    {
        // Check object on local mesh
        IOobject localIOobject
        (
            fieldName,
            localMesh.time().timeName(),
            cloud::prefix/cloudName,
            localMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (localIOobject.typeHeaderOk<IOField<Type>>(true))
        {
            IOField<Type> localField(localIOobject);

            const label offset = field.size();
            field.setSize(offset + localField.size());

            forAll(localField, j)
            {
                field[offset + j] = localField[j];
            }
        }
    }

    return tfield;
}


template<class Type>
Foam::tmp<Foam::CompactIOField<Foam::Field<Type>, Type>>
Foam::lagrangianReconstructor::reconstructFieldField
(
    const word& cloudName,
    const word& fieldName
)
{
    // Construct empty field on mesh
    auto tfield = tmp<CompactIOField<Field<Type>, Type>>::New
    (
        IOobject
        (
            fieldName,
            mesh_.time().timeName(),
            cloud::prefix/cloudName,
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        Field<Field<Type>>()
    );
    auto& field = tfield.ref();

    for (const fvMesh& localMesh : procMeshes_)
    {
        // Check object on local mesh
        IOobject localIOobject
        (
            fieldName,
            localMesh.time().timeName(),
            cloud::prefix/cloudName,
            localMesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if
        (
            localIOobject.typeHeaderOk<CompactIOField<Field<Type>, Type>>
            (
                false
            )
         || localIOobject.typeHeaderOk<IOField<Field<Type>>>(false)
        )
        {
            CompactIOField<Field<Type>, Type> localField(localIOobject);

            const label offset = field.size();
            field.setSize(offset + localField.size());

            forAll(localField, j)
            {
                field[offset + j] = localField[j];
            }
        }
    }

    return tfield;
}


template<class Type>
Foam::label Foam::lagrangianReconstructor::reconstructFields
(
    const word& cloudName,
    const UPtrList<const IOobject>& fieldObjects
)
{
    typedef IOField<Type> fieldType;

    label nFields = 0;

    for (const IOobject& io : fieldObjects)
    {
        if (io.isHeaderClass<fieldType>())
        {
            if (verbose_)
            {
                if (!nFields)
                {
                    Info<< "    Reconstructing lagrangian "
                        << fieldType::typeName << "s\n" << nl;
                }
                Info<< "        " << io.name() << endl;
            }
            ++nFields;

            reconstructField<Type>(cloudName, io.name())().write();
        }
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


template<class Type>
Foam::label Foam::lagrangianReconstructor::reconstructFields
(
    const word& cloudName,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef IOField<Type> fieldType;

    return reconstructFields<Type>
    (
        cloudName,
        (
            selectedFields.empty()
          ? objects.csorted<fieldType>()
          : objects.csorted<fieldType>(selectedFields)
        )
    );
}


template<class Type>
Foam::label Foam::lagrangianReconstructor::reconstructFieldFields
(
    const word& cloudName,
    const IOobjectList& objects,
    const wordRes& selectedFields
)
{
    typedef CompactIOField<Field<Type>, Type> fieldType;
    typedef IOField<Field<Type>> fieldTypeB;

    UPtrList<const IOobject> fieldObjects;

    if (selectedFields.empty())
    {
        fieldObjects.push_back(objects.csorted<fieldType>());
        fieldObjects.push_back(objects.csorted<fieldTypeB>());
    }
    else
    {
        fieldObjects.push_back(objects.csorted<fieldType>(selectedFields));
        fieldObjects.push_back(objects.csorted<fieldTypeB>(selectedFields));
    }

    Foam::sort(fieldObjects, nameOp<IOobject>());

    label nFields = 0;
    for (const IOobject& io : fieldObjects)
    {
        if (verbose_)
        {
            if (!nFields)
            {
                Info<< "    Reconstructing lagrangian "
                    << fieldType::typeName << "s\n" << nl;
            }
            Info<< "        " << io.name() << endl;
        }
        ++nFields;

        reconstructFieldField<Type>(cloudName, io.name())().write();
    }

    if (verbose_ && nFields) Info<< endl;
    return nFields;
}


// ************************************************************************* //
