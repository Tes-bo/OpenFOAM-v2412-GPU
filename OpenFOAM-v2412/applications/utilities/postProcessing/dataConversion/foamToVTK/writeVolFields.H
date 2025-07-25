/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Code chunk for reading volume fields from disk
    and write with vtk::internalWriter and vtk::patchWriter

\*---------------------------------------------------------------------------*/

#ifndef FoamToVTK_writeVolFields_H
#define FoamToVTK_writeVolFields_H

#include "readFields.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoField>
label writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar,
    objectRegistry* cache
)
{
    // Sanity test
    if (!internalWriter && patchWriters.empty()) return 0;

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        tmp<GeoField> tfield =
            getField<GeoField>(proxy, objects, fieldName, syncPar, cache);

        if (tfield)
        {
            // Internal
            if (internalWriter)
            {
                internalWriter->write(tfield());
            }

            // Boundary
            for (vtk::patchWriter& writer : patchWriters)
            {
                writer.write(tfield());
            }

            tfield.clear();
            ++count;
        }
    }

    return count;
}


template<class GeoField>
label writeVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar,
    objectRegistry* cache
)
{
    // Sanity test
    if
    (
        (!internalWriter || !pInterp)
     && (patchWriters.empty() || patchInterps.empty())
    )
    {
        return 0;
    }

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        tmp<GeoField> tfield =
            getField<GeoField>(proxy, objects, fieldName, syncPar, cache);

        if (tfield)
        {
            // Internal
            if (internalWriter && pInterp)
            {
                internalWriter->write(tfield(), *pInterp);
            }

            // Boundary
            label writeri = 0;
            for (vtk::patchWriter& writer : patchWriters)
            {
                if (patchInterps.test(writeri))
                {
                    writer.write(tfield(), patchInterps[writeri]);
                }
                ++writeri;
            }

            tfield.clear();
            ++count;
        }
    }

    return count;
}


label writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar,
    objectRegistry* cache
)
{
    // Sanity test
    if (!internalWriter && patchWriters.empty()) return 0;

    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter,                     \
            patchWriters,                       \
            proxy,                              \
            objects,                            \
            syncPar,                            \
            cache                               \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField);
    count += foamToVtk_WRITE_FIELD(volVectorField);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField);
    count += foamToVtk_WRITE_FIELD(volTensorField);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


label writeAllVolFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    const autoPtr<volPointInterpolation>& pInterp,

    UPtrList<vtk::patchWriter>& patchWriters,
    const UPtrList<PrimitivePatchInterpolation<primitivePatch>>& patchInterps,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar,
    objectRegistry* cache
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeVolFields<FieldType>               \
        (                                       \
            internalWriter, pInterp,            \
            patchWriters,   patchInterps,       \
            proxy,                              \
            objects,                            \
            syncPar,                            \
            cache                               \
        )


    label count = 0;
    count += foamToVtk_WRITE_FIELD(volScalarField);
    count += foamToVtk_WRITE_FIELD(volVectorField);
    count += foamToVtk_WRITE_FIELD(volSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(volSymmTensorField);
    count += foamToVtk_WRITE_FIELD(volTensorField);

    #undef foamToVtk_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
