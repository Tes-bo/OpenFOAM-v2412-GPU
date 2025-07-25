/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "faMesh.H"
#include "mapPolyMesh.H"
#include "MapFaFields.H"
#include "faMeshMapper.H"
#include "areaFields.H"
#include "edgeFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMesh::updateMesh(const mapPolyMesh& mpm)
{
    // if (!mpm.morphing())
    // {
    //     // No topo change
    //     return false;
    // }

    // Create fa mesh mapper, using the old mesh
    const faMeshMapper mapper(*this, mpm);

    // Rebuild mesh
    // ~~~~~~~~~~~~

    // Clear existing mesh data
    clearOut();

    // Set new labels
    faceLabels_ = mapper.areaMap().newFaceLabels();

    const uindirectPrimitivePatch& bp = patch();

    // Collect patch data
    const label nTotalEdges = bp.nEdges();
    const label nInternalEdges = bp.nInternalEdges();
    const label nBoundaryEdges = bp.nBoundaryEdges();
    const labelListList& edgeFaces = bp.edgeFaces();

    labelListList patchEdges(boundary_.size());

    // Special handling required for faces that have more than one edge
    // Each patch will be visited separately

    labelList edgeToPatch(nBoundaryEdges, -1);
    const labelList& newFaceLabelsMap = mapper.areaMap().newFaceLabelsMap();

    const labelListList& oldPatchEdgeFaces = mapper.oldPatchEdgeFaces();

    forAll(oldPatchEdgeFaces, patchI)
    {
        labelList& curPatchEdges = patchEdges[patchI];
        curPatchEdges.resize(nBoundaryEdges);
        label nCurPatchEdges = 0;

        // Note: it is possible to pick up the old-to-new boundary patch
        // mapping, but currently this is not done.  HJ, 13/Aug/2011

        // Make a fast lookup
        labelHashSet oldFaceLookup(oldPatchEdgeFaces[patchI]);

        for  (label edgeI = nInternalEdges; edgeI < nTotalEdges; ++edgeI)
        {
            if (edgeToPatch[edgeI - nInternalEdges] > -1)
            {
                // Edge already found; continue with the next one
                continue;
            }

            // Boundary edges will only have one face next to them
            const label oldFaceIndex = newFaceLabelsMap[edgeFaces[edgeI][0]];

            if (oldFaceIndex > -1)
            {
                // Old face exists.  See if it has got an edge in this patch
                if (oldFaceLookup.found(oldFaceIndex))
                {
                    // Face found, add it to the patch
                    curPatchEdges[nCurPatchEdges] = edgeI;
                    nCurPatchEdges++;

                    edgeToPatch[edgeI - nInternalEdges] = patchI;
                }
            }
        }

        // Collected all faces for the current patch
        curPatchEdges.setSize(nCurPatchEdges);
    }

    // Set new edges for all patches
    forAll(boundary_, patchI)
    {
        boundary_[patchI].resetEdges(patchEdges[patchI]);
    }

    setPrimitiveMeshData();

    // Create global mesh data
    if (Pstream::parRun())
    {
        globalData();
    }

    // Calculate topology for the patches (processor-processor comms etc.)
    boundary_.updateMesh();

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();


    // Map fields
    mapFields(mapper);

    // Map old areas
    mapOldAreas(mapper);

    // Update edge interpolation
    edgeInterpolation::movePoints();

    return;
}


void Foam::faMesh::mapFields(const faMeshMapper& mapper) const
{
    // Map areaFields in the objectRegistry
    MapGeometricFields<scalar, faPatchField, faMeshMapper, areaMesh>(mapper);
    MapGeometricFields<vector, faPatchField, faMeshMapper, areaMesh>(mapper);
    MapGeometricFields<sphericalTensor, faPatchField, faMeshMapper, areaMesh>
        (mapper);
    MapGeometricFields<symmTensor, faPatchField, faMeshMapper, areaMesh>
        (mapper);
    MapGeometricFields<tensor, faPatchField, faMeshMapper, areaMesh>(mapper);

    // Map edgeFields in the objectRegistry
    MapGeometricFields<scalar, faePatchField, faMeshMapper, edgeMesh>(mapper);
    MapGeometricFields<vector, faePatchField, faMeshMapper, edgeMesh>(mapper);
    MapGeometricFields<sphericalTensor, faePatchField, faMeshMapper, edgeMesh>
        (mapper);
    MapGeometricFields<symmTensor, faePatchField, faMeshMapper, edgeMesh>
        (mapper);
    MapGeometricFields<tensor, faePatchField, faMeshMapper, edgeMesh>(mapper);
}


void Foam::faMesh::mapOldAreas(const faMeshMapper& mapper) const
{
    if (S0Ptr_)
    {
        DebugInFunction << "Mapping old face areas." << endl;

        scalarField& S0 = *S0Ptr_;

        scalarField savedS0(S0);
        S0.setSize(nFaces());

        const labelList& faceMap = mapper.areaMap().newFaceLabelsMap();

        // Map existing old areas; for new faces set area to zero
        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] > -1)
            {
                S0[faceI] = savedS0[faceMap[faceI]];
            }
            else
            {
                S0[faceI] = 0;
            }
        }
    }

    if (S00Ptr_)
    {
        DebugInFunction << "Mapping old-old face areas." << endl;

        scalarField& S00 = *S00Ptr_;

        scalarField savedS00(S00);
        S00.setSize(nFaces());

        const labelList& faceMap = mapper.areaMap().newFaceLabelsMap();

        // Map old areas for existing faces; for new faces, set area to zero
        forAll(faceMap, faceI)
        {
            if (faceMap[faceI] > -1)
            {
                S00[faceI] = savedS00[faceMap[faceI]];
            }
            else
            {
                S00[faceI] = 0;
            }
        }
    }
}


// ************************************************************************* //
