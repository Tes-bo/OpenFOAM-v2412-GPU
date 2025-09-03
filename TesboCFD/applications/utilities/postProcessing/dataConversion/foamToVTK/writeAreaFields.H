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
    Code chunk for reading finite-area fields from disk
    and write with vtk::uindirectPatchGeoFieldsWriter

\*---------------------------------------------------------------------------*/

#ifndef FoamToVTK_writeAreaFields_H
#define FoamToVTK_writeAreaFields_H

#include "readFields.H"
#include "foamVtkIndPatchGeoFieldsWriter.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Writer type for finite-area mesh + fields
typedef vtk::uindirectPatchGeoFieldsWriter vtkWriterType_areaMesh;

template<class GeoField>
label writeAreaFields
(
    vtkWriterType_areaMesh& writer,
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    label count = 0;

    for (const word& fieldName : objects.sortedNames<GeoField>())
    {
        tmp<GeoField> tfield =
            getField<GeoField>(mesh, objects, fieldName, syncPar);

        if (tfield)
        {
            writer.write(tfield());

            tfield.clear();
            ++count;
        }
    }

    return count;
}


label writeAllAreaFields
(
    vtkWriterType_areaMesh& writer,
    const faMesh& mesh,
    const IOobjectList& objects,
    const bool syncPar
)
{
    #undef  foamToVtk_WRITE_FIELD
    #define foamToVtk_WRITE_FIELD(FieldType)    \
        writeAreaFields<FieldType>              \
        (                                       \
            writer,                             \
            mesh,                               \
            objects,                            \
            syncPar                             \
        )

    label count = 0;
    count += foamToVtk_WRITE_FIELD(areaScalarField);
    count += foamToVtk_WRITE_FIELD(areaVectorField);
    count += foamToVtk_WRITE_FIELD(areaSphericalTensorField);
    count += foamToVtk_WRITE_FIELD(areaSymmTensorField);
    count += foamToVtk_WRITE_FIELD(areaTensorField);

    #undef foamToVtk_WRITE_FIELD
    return count;
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
