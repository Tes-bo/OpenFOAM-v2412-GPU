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
    Code chunk for reading point fields from disk
    and write with vtk::internalWriter and vtk::patchWriter

\*---------------------------------------------------------------------------*/

#ifndef FoamToVTK_writePointFields_H
#define FoamToVTK_writePointFields_H

#include "readFields.H"
#include "foamVtkInternalWriter.H"
#include "foamVtkPatchWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class GeoField>
bool writePointField
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const tmp<GeoField>& tfield,
    const fvMeshSubsetProxy& proxy
)
{
    // Sanity test
    if (!internalWriter && patchWriters.empty()) return false;

    if (!tfield)
    {
        return false;
    }

    tmp<GeoField> tproxied;
    if (proxy.useSubMesh())
    {
        tproxied = proxy.interpolate(tfield());
        tfield.clear();
    }
    else
    {
        tproxied = tfield;
    }

    if (!tproxied)
    {
        // Or Error?
        return false;
    }


    const auto& field = tproxied();

    // Internal
    if (internalWriter)
    {
        internalWriter->write(field);
    }

    // Boundary
    for (vtk::patchWriter& writer : patchWriters)
    {
        writer.write(field);
    }


    tproxied.clear();

    return true;
}


template<class GeoField>
label writePointFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const typename GeoField::Mesh& ptMesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    // Sanity test
    if (!internalWriter && patchWriters.empty()) return 0;

    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        if
        (
            writePointField<GeoField>
            (
                internalWriter,
                patchWriters,
                getField<GeoField>(ptMesh, objects, fieldName, syncPar),
                proxy
            )
        )
        {
            ++count;
        }
    }

    return count;
}


label writeAllPointFields
(
    autoPtr<vtk::internalWriter>& internalWriter,
    UPtrList<vtk::patchWriter>& patchWriters,

    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const bool syncPar
)
{
    // Sanity test
    if (!internalWriter && patchWriters.empty()) return 0;

    const pointMesh& ptMesh = pointMesh::New(proxy.baseMesh());

    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writePointFields<FieldType>             \
        (                                       \
            internalWriter,                     \
            patchWriters,                       \
            proxy, ptMesh,                      \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(pointScalarField);
    count += foamToVtk_WRITE_FIELD(pointVectorField);
    count += foamToVtk_WRITE_FIELD(pointSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(pointSymmTensorField);
    count += foamToVtk_WRITE_FIELD(pointTensorField);

    #undef foamToVTK_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
