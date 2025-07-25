/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

InNamespace
    Foam

Description
    Read point fields from disk and write as ensight data

\*---------------------------------------------------------------------------*/

#ifndef ensight_writePointFields_H
#define ensight_writePointFields_H

#include "readFields.H"
#include "ensightMesh.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
bool writePointField
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const tmp<PointField<Type>>& tfield
)
{
    if (!tfield)
    {
        return false;
    }

    const auto& field = tfield();

    // PointData = true
    autoPtr<ensightFile> os = ensCase.newData<Type>(field.name(), true);

    bool wrote = ensightOutput::writePointField<Type>
    (
        os.ref(),
        field,
        ensMesh
    );

    tfield.clear();
    return wrote;
}


template<class Type>
label writePointFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    typedef PointField<Type> FieldType;

    const pointMesh& pMesh = pointMesh::New(ensMesh.mesh());

    label count = 0;

    for (const IOobject& io : objects.csorted<FieldType>())
    {
        if
        (
            writePointField<Type>
            (
                ensCase,
                ensMesh,
                getField<FieldType>(io, pMesh)
            )
        )
        {
            Info<< ' ' << io.name();
            ++count;
        }
    }

    return count;
}


label writeAllPointFields
(
    ensightCase& ensCase,
    const ensightMesh& ensMesh,
    const IOobjectList& objects
)
{
    #undef  ensight_WRITE_FIELD
    #define ensight_WRITE_FIELD(PrimitiveType)          \
        writePointFields<PrimitiveType>                 \
        (                                               \
            ensCase,                                    \
            ensMesh,                                    \
            objects                                     \
        )

    label count = 0;
    count += ensight_WRITE_FIELD(scalar);
    count += ensight_WRITE_FIELD(vector);
    count += ensight_WRITE_FIELD(sphericalTensor);
    count += ensight_WRITE_FIELD(symmTensor);
    count += ensight_WRITE_FIELD(tensor);

    #undef ensight_WRITE_FIELD
    return count;
}

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
