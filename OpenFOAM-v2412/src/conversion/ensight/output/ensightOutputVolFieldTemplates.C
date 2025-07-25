/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "ensightOutputVolField.H"
#include "ensightMesh.H"

#include "fvMesh.H"
#include "linear.H"
#include "volPointInterpolation.H"
#include "interpolation.H"
#include "processorFvPatch.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * *  Detail * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh
)
{
    bool parallel = Pstream::parRun();

    const fvMesh& mesh = vf.mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();

    // Write internalMesh and cellZones - sorted by index
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        ensightOutput::writeField(scratch, os, vf, part, parallel);
    }


    // Write patches - sorted by index
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        if (patchId < 0 || patchId >= pbm.size())
        {
            // Future handling of combined patches?
            continue;
        }

        const label patchStart = pbm[patchId].start();

        // Either use a flat boundary field for all patches,
        // or patch-local face ids

        // Operate on a copy
        ensightFaces localPart(part);

        // Change from global faceIds to patch-local
        localPart.decrFaceIds(patchStart);

        ensightOutput::writeField
        (
            scratch,
            os,
            vf.boundaryField()[patchId],
            localPart,
            parallel
        );
    }


    // No face zones data
    if (faceZoneParts.empty())
    {
        return true;
    }


    // Flat boundary field
    // similar to volPointInterpolation::flatBoundaryField()

    Field<Type> flat(pbm.nFaces(), Foam::zero{});

    forAll(vf.boundaryField(), patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const auto& pfld = vf.boundaryField()[patchi];

        // Note: restrict transcribing to actual size of the patch field
        // - handles "empty" patch type etc.

        SubList<Type> slice(flat, pfld.size(), pp.offset());

        if (isA<processorPolyPatch>(pp))
        {
            // Use average value for processor faces
            // own cell value = patchInternalField
            // nei cell value = evaluated boundary values
            slice = (0.5 * (pfld.patchInternalField() + pfld));
        }
        else if (!isA<emptyPolyPatch>(pp))
        {
            slice = pfld;
        }
    }


    // Interpolate cell values to faces
    // - only needed when exporting faceZones...

    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tsfld
        = Foam::linearInterpolate(vf);

    const auto& sfld = tsfld();


    // Local output buffer
    label maxLen = 0;

    forAllConstIters(faceZoneParts, iter)
    {
        maxLen = max(maxLen, iter.val().size());
    }

    DynamicField<Type> values(maxLen);


    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts.sortedToc())
    {
        const ensightFaces& part = faceZoneParts[zoneId];

        // Loop over face ids to store the needed field values
        // - internal faces use linear interpolation
        // - boundary faces use the corresponding patch value

        // Local copy of the field
        values.resize_nocopy(part.size());
        values = Zero;

        auto valIter = values.begin();

        for (const label faceId : part.faceIds())
        {
            *valIter =
            (
                mesh.isInternalFace(faceId)
              ? sfld[faceId]
              : flat[faceId - mesh.nInternalFaces()]
            );

            ++valIter;
        }

        // The field is already in the proper element order
        // - just need its corresponding sub-fields
        ensightOutput::Detail::writeFaceSubField
        (
            scratch, os, values, part, parallel
        );
    }

    return true;
}


template<class Type>
bool Foam::ensightOutput::writePointField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const GeometricField<Type, pointPatchField, pointMesh>& pf,
    const ensightMesh& ensMesh
)
{
    bool parallel = Pstream::parRun();

    const polyMesh& mesh = ensMesh.mesh();

    const Map<ensightCells>& cellZoneParts = ensMesh.cellZoneParts();
    const Map<ensightFaces>& faceZoneParts = ensMesh.faceZoneParts();
    const Map<ensightFaces>& boundaryParts = ensMesh.boundaryParts();

    //
    // Write internalMesh and cellZones - sorted by index
    //
    for (const label zoneId : cellZoneParts.sortedToc())
    {
        const ensightCells& part = cellZoneParts[zoneId];

        labelList uniquePointLabels;
        part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        // Skip if empty
        if
        (
            parallel
          ? returnReduceAnd(uniquePointLabels.empty())
          : uniquePointLabels.empty()
        )
        {
            continue;
        }

        ensightOutput::Detail::writeFieldComponents
        (
            scratch,
            os,
            ensightFile::coordinates,
            UIndirectList<Type>(pf.internalField(), uniquePointLabels),
            parallel
        );
    }


    //
    // Write patches - sorted by index
    //
    for (const label patchId : boundaryParts.sortedToc())
    {
        const ensightFaces& part = boundaryParts[patchId];

        labelList uniquePointLabels;
        part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        // Skip if empty
        if
        (
            parallel
          ? returnReduceAnd(uniquePointLabels.empty())
          : uniquePointLabels.empty()
        )
        {
            continue;
        }

        const auto& bfld = pf.boundaryField()[patchId];

        // Only valuePointPatchField is actually derived from Field
        const auto* vpp = isA<Field<Type>>(bfld);
        if (vpp)
        {
            // Require patch local indices, not mesh point labels.
            // But need to use polyPatch meshPointMap() to recover the
            // local indices since the ensight output will have jumbled
            // the face output order

            const polyPatch& pp = mesh.boundaryMesh()[patchId];

            for (label& pointi : uniquePointLabels)
            {
                pointi = pp.meshPointMap()[pointi];
            }

            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightFile::coordinates,
                UIndirectList<Type>(*vpp, uniquePointLabels),
                parallel
            );
        }
        else
        {
            ensightOutput::Detail::writeFieldComponents
            (
                scratch,
                os,
                ensightFile::coordinates,
                UIndirectList<Type>(pf.internalField(), uniquePointLabels),
                parallel
            );
        }
    }

    //
    // Write requested faceZones - sorted by index
    //
    for (const label zoneId : faceZoneParts.sortedToc())
    {
        const ensightFaces& part = faceZoneParts[zoneId];

        labelList uniquePointLabels;
        part.uniqueMeshPoints(mesh, uniquePointLabels, parallel);

        if (Pstream::master())
        {
            os.beginPart(part.index());
        }

        // Skip if empty
        if
        (
            parallel
          ? returnReduceAnd(uniquePointLabels.empty())
          : uniquePointLabels.empty()
        )
        {
            continue;
        }

        // CAVEAT - does not properly handle valuePointPatchField,
        // uses internalField only

        ensightOutput::Detail::writeFieldComponents
        (
            scratch,
            os,
            ensightFile::coordinates,
            UIndirectList<Type>(pf.internalField(), uniquePointLabels),
            parallel
        );
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::ensightOutput::writeVolField
(
    ensightOutput::floatBufferType& scratch,
    ensightFile& os,
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const ensightMesh& ensMesh,
    const bool nodeValues
)
{
    if (nodeValues)
    {
        tmp<GeometricField<Type, pointPatchField, pointMesh>> pfld
        (
            volPointInterpolation::New(vf.mesh()).interpolate(vf)
        );
        pfld.ref().checkOut();
        pfld.ref().rename(vf.name());

        return ensightOutput::writePointField<Type>(scratch, os, pfld, ensMesh);
    }

    return ensightOutput::writeVolField<Type>(scratch, os, vf, ensMesh);
}


// ************************************************************************* //
