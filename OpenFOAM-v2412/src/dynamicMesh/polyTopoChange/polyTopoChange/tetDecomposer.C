/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020,2022,2024 OpenCFD Ltd.
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

#include "tetDecomposer.H"
#include "meshTools.H"
#include "polyMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "OFstream.H"
#include "edgeHashes.H"
#include "syncTools.H"
#include "dummyTransform.H"
#include "triangle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(tetDecomposer, 0);
}

const Foam::Enum
<
    Foam::tetDecomposer::decompositionType
>
Foam::tetDecomposer::decompositionTypeNames
({
    { decompositionType::FACE_CENTRE_TRIS,  "faceCentre" },
    { decompositionType::FACE_DIAG_TRIS, "faceDiagonal" },
    { decompositionType::PYRAMID, "pyramid" },
    { decompositionType::FACE_DIAG_QUADS, "faceDiagonalQuads" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::tetDecomposer::modifyFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const label patchi,
    const label zoneI,
    const bool zoneFlip
) const
{
    if (own == nei)
    {
        FatalErrorInFunction
            << "Problem own:" << own
            << " nei:" << nei << exit(FatalError);
    }

    // First usage of face. Modify.
    if (nei == -1 || own < nei)
    {
        meshMod.modifyFace
        (
            f,                          // modified face
            facei,                      // label of face
            own,                        // owner
            nei,                        // neighbour
            false,                      // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        meshMod.modifyFace
        (
            f.reverseFace(),            // modified face
            facei,                      // label of face
            nei,                        // owner
            own,                        // neighbour
            true,                       // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


void Foam::tetDecomposer::addFace
(
    polyTopoChange& meshMod,
    const face& f,
    const label facei,
    const label own,
    const label nei,
    const label masterPointID,
    const label masterEdgeID,
    const label masterFaceID,
    const label patchi,
    const label zoneI,
    const bool zoneFlip
) const
{
    if (own == nei)
    {
        FatalErrorInFunction
            << "Problem own:" << own
            << " nei:" << nei << exit(FatalError);
    }

    // Second or more usage of face. Add.
    if (nei == -1 || own < nei)
    {
        //const label newFacei =
        meshMod.addFace
        (
            f,                          // modified face
            own,                        // owner
            nei,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            false,                      // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            zoneFlip                    // face flip in zone
        );
    }
    else
    {
        //const label newFacei =
        meshMod.addFace
        (
            f.reverseFace(),            // modified face
            nei,                        // owner
            own,                        // neighbour
            masterPointID,              // master point
            masterEdgeID,               // master edge
            masterFaceID,               // master face
            true,                       // face flip
            patchi,                     // patch for face
            zoneI,                      // zone for face
            !zoneFlip                   // face flip in zone
        );
    }
}


// Work out triangle index given the starting vertex in the face
Foam::label Foam::tetDecomposer::triIndex(const label facei, const label fp)
const
{
    const face& f = mesh_.faces()[facei];
    const label fp0 = max(0, mesh_.tetBasePtIs()[facei]);

    // Work out triangle index on this face. Assume 'fan' triangulation starting
    // from fp0. E.g. 6 vertices on face -> 4 triangles. First and last triangle
    // use consecutive vertices
    //
    //    fp    | triangle | vertices
    //    ------|----------|---------
    //    fp0   |  0       | fp0,fp0+1,fp0+2
    //    fp0+1 |  0       |  ,,
    //    fp0+2 |  1       | fp0,fp0+2,fp0+3
    //    fp0+3 |  2       | fp0,fp0+3,fp0+4
    //    fp0+4 |  3       | fp0,fp0+4,fp0+3
    //    fp0+5 |  3       |  ,,


    // Work out triangle index on this face
    label thisTrii;
    if (fp == fp0)
    {
        thisTrii = 0;
    }
    else if (fp == f.fcIndex(fp0))
    {
        thisTrii = 0;
    }
    else if (fp == f.rcIndex(fp0))
    {
        thisTrii = f.size()-3;
    }
    else if (fp < fp0)
    {
        const label fpB = fp+f.size();
        thisTrii = (fpB-fp0-1);
    }
    else
    {
        thisTrii = (fp-fp0-1);
    }
    return thisTrii;
}


void Foam::tetDecomposer::splitBoundaryFaces
(
    List<faceList>& boundaryQuads,
    List<faceList>& boundaryTris
) const
{
    // Work space
    faceList quadFaces(1000);
    faceList triFaces(1000);

    const auto& pbm = mesh_.boundaryMesh();
    for (const auto& pp : pbm)
    {
        if (pp.coupled() && refCast<const coupledPolyPatch>(pp).owner())
        {
            forAll(pp, i)
            {
                const label facei = pp.start()+i;
                const face& meshFace = pp[i];

                if (meshFace.size() > 4)
                {
                    label trii = 0;
                    label quadi = 0;
                    meshFace.trianglesQuads
                    (
                        mesh_.points(),
                        trii,
                        quadi,
                        triFaces,
                        quadFaces
                    );

                    const label bFacei = facei-mesh_.nInternalFaces();

                    {
                        auto& faces = boundaryTris[bFacei];
                        faces.setSize(trii);
                        // Express as relative w.r.t. 0th vertex
                        for (label i = 0; i < trii; i++)
                        {
                            const auto& f = triFaces[i];
                            auto& verts = faces[i];
                            verts.setSize(f.size());
                            forAll(f, fp)
                            {
                                verts[fp] = meshFace.find(f[fp]);
                            }
                        }
                    }
                    {
                        auto& faces = boundaryQuads[bFacei];
                        faces.setSize(quadi);
                        // Express as relative w.r.t. 0th vertex
                        for (label i = 0; i < quadi; i++)
                        {
                            const auto& f = quadFaces[i];
                            auto& verts = faces[i];
                            verts.setSize(f.size());
                            forAll(f, fp)
                            {
                                verts[fp] = meshFace.find(f[fp]);
                            }
                        }
                    }
                }
            }
        }
    }

    // Send coupled side indices to neighbour. Note: could also do re-indexing
    // here...
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        boundaryTris,
        [](faceList& result, const faceList& input)
        {
            if (!result.size())
            {
                result = input;
            }
        },
        dummyTransform()
    );
    syncTools::syncBoundaryFaceList
    (
        mesh_,
        boundaryQuads,
        [](faceList& result, const faceList& input)
        {
            if (!result.size())
            {
                result = input;
            }
        },
        dummyTransform()
    );
}


void Foam::tetDecomposer::relativeIndicesToFace
(
    const bool doFlip,
    const face& meshFace,
    const faceList& indexLists,
    faceList& faces
) const
{
    //faces.setSize(indexLists.size());

    if (!doFlip)
    {
        forAll(indexLists, facei)
        {
            const auto& verts = indexLists[facei];
            auto& f = faces[facei];
            f.setSize(verts.size());

            forAll(verts, fp)
            {
                f[fp] = meshFace[verts[fp]];
            }
        }
    }
    else
    {
        forAllReverse(indexLists, facei)
        {
            const auto& verts = indexLists[facei];
            auto& f = faces[facei];
            f.setSize(verts.size());

            // - 0th vertex matches; walking order opposite
            // - assemble in opposite order so as to flip normal
            label destFp = verts.size()-1;
            forAll(verts, fp)
            {
                if (verts[fp] == 0)
                {
                    f[destFp] = meshFace[0];
                }
                else
                {
                    f[destFp] = meshFace[meshFace.size()-verts[fp]];
                }
                --destFp;
            }
        }
    }
}


void Foam::tetDecomposer::splitFace
(
    const List<faceList>& boundaryQuads,
    const List<faceList>& boundaryTris,
    const label facei,
    const label patchi,
    label& quadi,
    faceList& quadFaces,
    label& trii,
    faceList& triFaces
) const
{
    // Split face into quads (in quadFaces) and tris (in triFaces)

    const face& f = mesh_.faces()[facei];

    trii = 0;
    quadi = 0;
    if (patchi != -1 && mesh_.boundaryMesh()[patchi].coupled())
    {
        const auto& pp = mesh_.boundaryMesh()[patchi];
        const bool owner =
            refCast<const coupledPolyPatch>(pp).owner();
        const label bFacei = facei-mesh_.nInternalFaces();

        // Triangles
        trii = boundaryTris[bFacei].size();
        relativeIndicesToFace
        (
            !owner,
            f,
            boundaryTris[bFacei],
            triFaces
        );

        // Quads
        quadi = boundaryQuads[bFacei].size();
        relativeIndicesToFace
        (
            !owner,
            f,
            boundaryQuads[bFacei],
            quadFaces
        );
    }
    else if (f.size() == 4)
    {
        quadFaces[quadi++] = f;
    }
    else
    {
        f.trianglesQuads
        (
            mesh_.points(),
            trii,
            quadi,
            triFaces,
            quadFaces
        );
    }
}


void Foam::tetDecomposer::splitFace
(
    const List<faceList>& boundaryQuads,
    const List<faceList>& boundaryTris,
    const label facei,
    const label patchi,
    faceList& quadFaces,
    faceList& triFaces,
    faceList& subFaces
) const
{
    const face& f = mesh_.faces()[facei];

    if (f.size() <= 4)
    {
        subFaces.resize_nocopy(1);
        subFaces[0] = f;
    }
    else
    {
        label trii;
        label quadi;
        splitFace
        (
            boundaryQuads,
            boundaryTris,
            facei,
            patchi,
            quadi,
            quadFaces,
            trii,
            triFaces
        );

        // Collect into single face list
        subFaces.resize_nocopy(quadi+trii);
        {
            label subFacei = 0;
            for (label i = 0; i < quadi; i++)
            {
                subFaces[subFacei++] = quadFaces[i];
            }
            for (label i = 0; i < trii; i++)
            {
                subFaces[subFacei++] = triFaces[i];
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetDecomposer::tetDecomposer(const polyMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::tetDecomposer::setRefinement
(
    const decompositionType decomposeType,
    const bitSet& decomposeCell,
    polyTopoChange& meshMod
)
{
    // Determine for every face whether it borders a cell that is decomposed
    bitSet decomposeFace(mesh_.nFaces());
    {
        for (label facei = 0; facei < mesh_.nInternalFaces(); facei++)
        {
            label own = mesh_.faceOwner()[facei];
            label nei = mesh_.faceNeighbour()[facei];
            if (decomposeCell[own] || decomposeCell[nei])
            {
                decomposeFace[facei] = true;
            }
        }

        boolList neiDecomposeCell(mesh_.nBoundaryFaces());
        forAll(neiDecomposeCell, bFacei)
        {
            label facei = mesh_.nInternalFaces()+bFacei;
            label own = mesh_.faceOwner()[facei];
            neiDecomposeCell[bFacei] = decomposeCell[own];
        }
        syncTools::swapBoundaryFaceList(mesh_, neiDecomposeCell);

        for
        (
            label facei = mesh_.nInternalFaces();
            facei < mesh_.nFaces();
            facei++
        )
        {
            label own = mesh_.faceOwner()[facei];
            label bFacei = facei-mesh_.nInternalFaces();
            if (decomposeCell[own] || neiDecomposeCell[bFacei])
            {
                decomposeFace[facei] = true;
            }
        }
    }
    setRefinement(decomposeType, decomposeCell, decomposeFace, meshMod);
}


void Foam::tetDecomposer::setRefinement
(
    const decompositionType decomposeType,
    const bitSet& decomposeCell,
    const bitSet& decomposeFace,
    polyTopoChange& meshMod
)
{
    if (debug)
    {
        // Check that current mesh makes sense
        const auto& faces = mesh_.faces();
        labelList nVerts(mesh_.nFaces());
        forAll(faces, facei)
        {
            nVerts[facei] = faces[facei].size();
        }
        syncTools::swapFaceList(mesh_, nVerts);
        forAll(nVerts, facei)
        {
            if (nVerts[facei] != faces[facei].size())
            {
                FatalErrorInFunction<< "problem with nVerts"
                    << exit(FatalError);
            }
        }
    }
    if (debug)
    {
        // Check that decomposeFace makes sense
        bitSet newDecomposeFace(decomposeFace);
        syncTools::swapFaceList(mesh_, newDecomposeFace);
        forAll(newDecomposeFace, facei)
        {
            if (decomposeFace[facei] != newDecomposeFace[facei])
            {
                FatalErrorInFunction<< "problem with decomposeFace"
                    << exit(FatalError);
            }
        }
    }


    // For diagonal splitting : send over diagonals since two sides of
    // coupled face might make different decisions
    // TBD: probably also for FACE_DIAG_TRIS ?
    faceList quadFaces;
    faceList triFaces;
    List<faceList> boundaryQuads;
    List<faceList> boundaryTris;
    faceList subFaces;

    if (decomposeType == FACE_DIAG_QUADS)
    {
        // Pre-calculate coupled faces triangulation. Store as indices w.r.t.
        // vertex 0.
        // Note: could triangulate all faces here but trying to save some
        // memory so only doing:
        //  - faces on proc boundaries
        //  - faces > 4 verts
        quadFaces.resize_nocopy(1000);
        triFaces.resize_nocopy(1000);
        boundaryQuads.resize_nocopy(mesh_.nBoundaryFaces());
        boundaryTris.resize_nocopy(mesh_.nBoundaryFaces());
        splitBoundaryFaces(boundaryQuads, boundaryTris);
    }


    cellToPoint_.setSize(mesh_.nCells(), -1);
    forAll(mesh_.cellCentres(), celli)
    {
        if (decomposeCell[celli])
        {
            // Any point on the cell
            label masterPointi = mesh_.faces()[mesh_.cells()[celli][0]][0];

            cellToPoint_[celli] = meshMod.addPoint
            (
                mesh_.cellCentres()[celli],
                masterPointi,
                -1,
                true
            );
        }
    }


    // Add face centre points
    if (decomposeType == FACE_CENTRE_TRIS)
    {
        faceToPoint_.setSize(mesh_.nFaces(), -1);
        forAll(mesh_.faceCentres(), facei)
        {
            if (decomposeFace[facei])
            {
                // Any point on the face
                const label masterPointi = mesh_.faces()[facei][0];

                faceToPoint_[facei] = meshMod.addPoint
                (
                    mesh_.faceCentres()[facei],
                    masterPointi,
                    -1,
                    true
                );
            }
        }
    }


    // Per face, per point (faceCentre) or triangle (faceDiag) the (existing
    // or added) cell on either side
    faceOwnerCells_.setSize(mesh_.nFaces());
    faceNeighbourCells_.setSize(mesh_.nFaces());

    if (decomposeType == FACE_CENTRE_TRIS)
    {
        forAll(faceOwnerCells_, facei)
        {
            if (decomposeFace[facei])
            {
                const face& f = mesh_.faces()[facei];
                faceOwnerCells_[facei].setSize(f.size(), -1);
                faceNeighbourCells_[facei].setSize(f.size(), -1);
            }
            else
            {
                faceOwnerCells_[facei].setSize(1, -1);
                faceNeighbourCells_[facei].setSize(1, -1);
            }
        }
    }
    else if (decomposeType == FACE_DIAG_TRIS)
    {
        // Force construction of diagonal decomposition
        (void)mesh_.tetBasePtIs();

        forAll(faceOwnerCells_, facei)
        {
            if (decomposeFace[facei])
            {
                const face& f = mesh_.faces()[facei];
                faceOwnerCells_[facei].setSize(f.size()-2, -1);
                faceNeighbourCells_[facei].setSize(f.size()-2, -1);
            }
            else
            {
                faceOwnerCells_[facei].setSize(1, -1);
                faceNeighbourCells_[facei].setSize(1, -1);
            }
        }
    }
    else if (decomposeType == FACE_DIAG_QUADS)
    {
        // Note: sizing according to the number of sub-faces, not according to
        //       f.size(). Plus: less storage. Min: much harder to look up
        //       cell corresponding to face+edge

        // Do coupled faces first - do not use local face triangulation
        // but use coupled version instead
        const auto& pbm = mesh_.boundaryMesh();
        for (const auto& pp : pbm)
        {
            if (pp.coupled())
            {
                forAll(pp, i)
                {
                    const label facei = pp.start()+i;
                    if (decomposeFace[facei])
                    {
                        const label bFacei = pp.offset()+i;
                        const label n =
                            boundaryQuads[bFacei].size()
                          + boundaryTris[bFacei].size();

                        faceOwnerCells_[facei].setSize(n, -1);
                        faceNeighbourCells_[facei].setSize(n, -1);
                    }
                }
            }
        }

        // Do all other faces
        forAll(faceOwnerCells_, facei)
        {
            if (faceOwnerCells_[facei].empty())
            {
                if (decomposeFace[facei])
                {
                    const face& f = mesh_.faces()[facei];

                    // Convention: quads first, followed by triangles

                    label nTris = 0;
                    label nQuads = 0;
                    const label nSubFaces = f.nTrianglesQuads
                    (
                        mesh_.points(),
                        nTris,
                        nQuads
                    );

                    faceOwnerCells_[facei].setSize(nSubFaces, -1);
                    faceNeighbourCells_[facei].setSize(nSubFaces, -1);
                }
                else
                {
                    faceOwnerCells_[facei].setSize(1, -1);
                    faceNeighbourCells_[facei].setSize(1, -1);
                }
            }
        }
    }
    else
    {
        forAll(faceOwnerCells_, facei)
        {
            faceOwnerCells_[facei].setSize(1, -1);
            faceNeighbourCells_[facei].setSize(1, -1);
        }
    }


    // Add internal cells. Note: done in same order as pyramid triangle
    // creation later to maintain same ordering.
    forAll(mesh_.cells(), celli)
    {
        const cell& cFaces = mesh_.cells()[celli];

        // Whether cell has already been modified (all other cells get added)
        bool modifiedCell = false;

        forAll(cFaces, cFacei)
        {
            label facei = cFaces[cFacei];

            // Get reference to either owner or neighbour
            labelList& added =
            (
                (mesh_.faceOwner()[facei] == celli)
              ? faceOwnerCells_[facei]
              : faceNeighbourCells_[facei]
            );

            if (decomposeCell[celli])
            {
                if (decomposeType == FACE_CENTRE_TRIS)
                {
                    forAll(added, i)
                    {
                        if (!modifiedCell)
                        {
                            // Reuse cell itself
                            added[i] = celli;
                            modifiedCell = true;
                        }
                        else
                        {
                            added[i] = meshMod.addCell
                            (
                                -1,     // masterPoint
                                -1,     // masterEdge
                                -1,     // masterFace
                                celli,  // masterCell
                                mesh_.cellZones().whichZone(celli)
                            );
                        }
                    }
                }
                else if
                (
                    decomposeType == FACE_DIAG_TRIS
                 || decomposeType == FACE_DIAG_QUADS
                )
                {
                    forAll(added, subi)
                    {
                        if (!modifiedCell)
                        {
                            // Reuse cell itself
                            added[subi] = celli;
                            modifiedCell = true;
                        }
                        else
                        {
                            added[subi] = meshMod.addCell
                            (
                                -1,     // masterPoint
                                -1,     // masterEdge
                                -1,     // masterFace
                                celli,  // masterCell
                                mesh_.cellZones().whichZone(celli)
                            );
                        }
                    }
                }
                else // if (decomposeType == PYRAMID)
                {
                    // Pyramidal decomposition.
                    // Assign same cell to all slots

                    if (added.size() != 1)
                    {
                        FatalErrorInFunction
                            << "For cell:" << celli
                            << " at:" << mesh_.cellCentres()[celli]
                            << " face:" << facei
                            << " at:" << mesh_.faceCentres()[facei]
                            << " added cells:" << added
                            << exit(FatalError);
                    }


                    if (!modifiedCell)
                    {
                        // Reuse cell itself
                        added = celli;
                        modifiedCell = true;
                    }
                    else
                    {
                        added = meshMod.addCell
                        (
                            -1,     // masterPoint
                            -1,     // masterEdge
                            -1,     // masterFace
                            celli,  // masterCell
                            mesh_.cellZones().whichZone(celli)
                        );
                    }
                }
            }
            else
            {
                // All vertices/triangles address to original cell
                added = celli;
            }
        }
    }



    // Add split faces
    face triangle(3);

    forAll(mesh_.faces(), facei)
    {
        label own = mesh_.faceOwner()[facei];
        const auto& addedOwn = faceOwnerCells_[facei];
        const auto& addedNei = faceNeighbourCells_[facei];
        const face& f = mesh_.faces()[facei];
        //const point& fc = mesh_.faceCentres()[facei];

        label nei = -1;
        label patchi = -1;
        if (facei >= mesh_.nInternalFaces())
        {
            patchi = mesh_.boundaryMesh().whichPatch(facei);
        }
        else
        {
            nei = mesh_.faceNeighbour()[facei];
        }

        label zonei = mesh_.faceZones().whichZone(facei);
        bool zoneFlip = false;
        if (zonei != -1)
        {
            const faceZone& fz = mesh_.faceZones()[zonei];
            zoneFlip = fz.flipMap()[fz.whichFace(facei)];
        }

        if (decomposeType == FACE_CENTRE_TRIS && decomposeFace[facei])
        {
            forAll(f, fp)
            {
                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    triangle[0] = f[fp];
                    triangle[1] = f[f.fcIndex(fp)];
                    triangle[2] = faceToPoint_[facei];

                    const label newOwn
                    (
                        addedOwn.size() == 1
                      ? addedOwn[0]
                      : addedOwn[fp]
                    );
                    const label newNei
                    (
                        addedNei.size() == 1
                      ? addedNei[0]
                      : addedNei[fp]
                    );

                    if (fp == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            -1,                 //point
                            -1,                 //edge
                            facei,              //face
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - to cell centre
                if (decomposeCell[own])
                {
                    label newOwn = addedOwn[f.rcIndex(fp)];
                    label newNei = addedOwn[fp];

                    triangle[0] = f[fp];
                    triangle[1] = cellToPoint_[own];
                    triangle[2] = faceToPoint_[facei];

                    addFace
                    (
                        meshMod,
                        triangle,
                        facei,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }
                // 2b. Within neighbour cell - to cell centre
                if (facei < mesh_.nInternalFaces() && decomposeCell[nei])
                {
                    label newOwn = addedNei[f.rcIndex(fp)];
                    label newNei = addedNei[fp];

                    triangle[0] = f[fp];
                    triangle[1] = faceToPoint_[facei];
                    triangle[2] = cellToPoint_[nei];

                    addFace
                    (
                        meshMod,
                        triangle,
                        facei,
                        newOwn,
                        newNei,
                        f[fp],      //point
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }
            }
        }
        else if (decomposeType == FACE_DIAG_TRIS && decomposeFace[facei])
        {
            label fp0 = max(mesh_.tetBasePtIs()[facei], 0);
            label fp = f.fcIndex(fp0);

            for (label trii = 0; trii < f.size()-2; trii++)
            {
                label nextTri = trii+1;
                if (nextTri >= f.size()-2)
                {
                    nextTri -= f.size()-2;
                }
                label nextFp = f.fcIndex(fp);


                // Triangle trii consisiting of f[fp0], f[fp], f[nextFp]


                // 1. Front triangle (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    const label newOwn
                    (
                        addedOwn.size() == 1
                      ? addedOwn[0]
                      : addedOwn[trii]
                    );
                    const label newNei
                    (
                        addedNei.size() == 1
                      ? addedNei[0]
                      : addedNei[trii]
                    );

                    triangle[0] = f[fp0];
                    triangle[1] = f[fp];
                    triangle[2] = f[nextFp];

                    if (trii == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            -1,                 //point
                            -1,                 //edge
                            facei,              //face
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                }


                // 2. Within owner cell - diagonal to cell centre
                if (trii < f.size()-3)
                {
                    if (decomposeCell[own])
                    {
                        label newOwn = addedOwn[trii];
                        label newNei = addedOwn[nextTri];

                        triangle[0] = f[fp0];
                        triangle[1] = f[nextFp];
                        triangle[2] = cellToPoint_[own];

                        addFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            f[fp],      //point
                            -1,         //edge
                            -1,         //face
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }

                    // 2b. Within neighbour cell - to cell centre
                    if (facei < mesh_.nInternalFaces() && decomposeCell[nei])
                    {
                        label newOwn = addedNei[trii];
                        label newNei = addedNei[nextTri];

                        triangle[0] = f[nextFp];
                        triangle[1] = f[fp0];
                        triangle[2] = cellToPoint_[nei];

                        addFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            newOwn,
                            newNei,
                            f[fp],      //point
                            -1,         //edge
                            -1,         //face
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }
                }


                fp = nextFp;
            }
        }
        else if (decomposeType == FACE_DIAG_QUADS && decomposeFace[facei])
        {
            // Decompose face into subFaces (quads followed by any triangles)
            splitFace
            (
                boundaryQuads,
                boundaryTris,
                facei,
                patchi,
                quadFaces,  // work space
                triFaces,   // work space
                subFaces
            );

            EdgeMap<labelPair> edgeFaces(subFaces.size());

            forAll(subFaces, subFacei)
            {
                const face& subF = subFaces[subFacei];

                // 1. Front quad (decomposition of face itself)
                //    (between owner and neighbour cell)
                {
                    const label newOwn
                    (
                        addedOwn.size() == 1
                      ? addedOwn[0]
                      : addedOwn[subFacei]
                    );
                    const label newNei
                    (
                        addedNei.size() == 1
                      ? addedNei[0]
                      : addedNei[subFacei]
                    );

                    if (subFacei == 0)
                    {
                        modifyFace
                        (
                            meshMod,
                            subF,
                            facei,
                            newOwn,
                            newNei,
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                    else
                    {
                        addFace
                        (
                            meshMod,
                            subF,
                            facei,
                            newOwn,
                            newNei,
                            -1,                 //point
                            -1,                 //edge
                            facei,              //face
                            patchi,
                            zonei,
                            zoneFlip
                        );
                    }
                }

                // Populate edge-faces (note: in order of increasing subFacei
                // and hence in order of added cells)
                forAll(subF, fp)
                {
                    const edge e(subF.edge(fp));
                    auto eFnd = edgeFaces.find(e);
                    if (!eFnd)
                    {
                        edgeFaces.insert(e, labelPair(subFacei, -1));
                    }
                    else
                    {
                        auto& eFaces = eFnd();
                        if (eFaces[1] != -1)
                        {
                            FatalErrorInFunction << "edgeFaces:" << edgeFaces
                            << exit(FatalError);
                        }
                        eFaces[1] = subFacei;
                    }
                }
            }

            // Get diagonals
            forAllConstIters(edgeFaces, iter)
            {
                const auto& e = iter.key();
                const auto& eFaces = iter();

                if (eFaces.find(-1) != -1)
                {
                    // Open edge
                    //Pout<< "for face:" << facei
                    //    << " ignoring open edge:" << e << endl;
                    continue;
                }

                if (decomposeCell[own])
                {
                    const label newOwn(addedOwn[eFaces[0]]);
                    const label newNei(addedOwn[eFaces[1]]);

                    if (newNei < newOwn) FatalErrorInFunction << "problem"
                        << exit(FatalError);

                    // Point out of owner side
                    triangle[0] = e[1];
                    triangle[1] = e[0];
                    triangle[2] = cellToPoint_[own];

                    addFace
                    (
                        meshMod,
                        triangle,
                        facei,
                        newOwn,
                        newNei,
                        e[0],       //point ? or edge ? or face ?
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }

                // 2b. Within neighbour cell - to cell centre
                if
                (
                    facei < mesh_.nInternalFaces()
                 && decomposeCell[nei]
                )
                {
                    const label newOwn(addedNei[eFaces[0]]);
                    const label newNei(addedNei[eFaces[1]]);

                    if (newNei < newOwn) FatalErrorInFunction << "problem"
                        << exit(FatalError);

                    triangle[0] = e[0];
                    triangle[1] = e[1];
                    triangle[2] = cellToPoint_[nei];

                    addFace
                    (
                        meshMod,
                        triangle,
                        facei,
                        newOwn,
                        newNei,
                        e[0],       //point ? or edge ? or face ?
                        -1,         //edge
                        -1,         //face
                        -1,         //patchi
                        -1,         //zone
                        false
                    );
                }
            }
        }
        else
        {
            // No decomposition. Use zero'th slot.
            modifyFace
            (
                meshMod,
                f,
                facei,
                addedOwn[0],
                addedNei[0],
                patchi,
                zonei,
                zoneFlip
            );
        }
    }


    // Add triangles to the cell centre for all edges to form the pyramids
    EdgeMap<label> edgeToFace;

    forAll(mesh_.cells(), celli)
    {
        if (decomposeCell[celli])
        {
            const cell& cFaces = mesh_.cells()[celli];

            edgeToFace.clear();

            forAll(cFaces, cFacei)
            {
                const label facei = cFaces[cFacei];
                const label patchi = mesh_.boundaryMesh().whichPatch(facei);
                const face& f = mesh_.faces()[facei];

                // Loop over edges of face. Avoid constructing faceEdges for
                // whole mesh.
                forAll(f, fp)
                {
                    label p0 = f[fp];
                    label p1 = f[f.fcIndex(fp)];
                    const edge e(p0, p1);

                    EdgeMap<label>::const_iterator edgeFnd = edgeToFace.find(e);
                    if (edgeFnd == edgeToFace.end())
                    {
                        edgeToFace.insert(e, facei);
                    }
                    else
                    {
                        // Found the other face on the edge.
                        label otherFacei = edgeFnd();
                        const face& otherF = mesh_.faces()[otherFacei];

                        // Found the other face on the edge. Note that since
                        // we are looping in the same order the tets added for
                        // otherFacei will be before those of facei

                        label otherFp = otherF.find(p0);
                        if (otherF.nextLabel(otherFp) == p1)
                        {
                            // ok. otherFp is first vertex of edge.
                        }
                        else if (otherF.prevLabel(otherFp) == p1)
                        {
                            otherFp = otherF.rcIndex(otherFp);
                        }
                        else
                        {
                            FatalErrorInFunction
                                << "problem." << abort(FatalError);
                        }


                        // Triangle from edge to cell centre
                        if (mesh_.faceOwner()[facei] == celli)
                        {
                            triangle[0] = p0;
                            triangle[1] = p1;
                            triangle[2] = cellToPoint_[celli];
                        }
                        else
                        {
                            triangle[0] = p1;
                            triangle[1] = p0;
                            triangle[2] = cellToPoint_[celli];
                        }

                        // Determine tets on either side

                        const auto& thisCells =
                        (
                            (mesh_.faceOwner()[facei] == celli)
                          ? faceOwnerCells_[facei]
                          : faceNeighbourCells_[facei]
                        );
                        const auto& otherCells =
                        (
                            (mesh_.faceOwner()[otherFacei] == celli)
                          ? faceOwnerCells_[otherFacei]
                          : faceNeighbourCells_[otherFacei]
                        );

                        label thisTet = -1;
                        label otherTet = -1;

                        if (decomposeType == FACE_CENTRE_TRIS)
                        {
                            if (thisCells.size() == 1)
                            {
                                thisTet = thisCells[0];
                            }
                            else
                            {
                                thisTet = thisCells[fp];
                            }

                            if (otherCells.size() == 1)
                            {
                                otherTet = otherCells[0];
                            }
                            else
                            {
                                otherTet = otherCells[otherFp];
                            }
                        }
                        else if (decomposeType == FACE_DIAG_TRIS)
                        {
                            if (thisCells.size() == 1)
                            {
                                thisTet = thisCells[0];
                            }
                            else
                            {
                                thisTet = thisCells[triIndex(facei, fp)];
                            }

                            if (otherCells.size() == 1)
                            {
                                otherTet = otherCells[0];
                            }
                            else
                            {
                                otherTet =
                                    otherCells[triIndex(otherFacei, otherFp)];
                            }
                        }
                        else if (decomposeType == FACE_DIAG_QUADS)
                        {
                            // Get added cell on facei
                            if (thisCells.size() == 1)
                            {
                                thisTet = thisCells[0];
                            }
                            else
                            {
                                // We have not stored the face-split so
                                // don't know yet the correspondence. Maybe
                                // store to avoid re-doing split?
                                splitFace
                                (
                                    boundaryQuads,
                                    boundaryTris,
                                    facei,
                                    patchi,
                                    quadFaces,  // work space
                                    triFaces,   // work space
                                    subFaces
                                );

                                // Find the edge. Should be only on one of the
                                // subfaces.
                                forAll(subFaces, subFacei)
                                {
                                    const auto& subF = subFaces[subFacei];
                                    if (subF.find(e) != -1)
                                    {
                                        thisTet = thisCells[subFacei];
                                        break;
                                    }
                                }
                                if (thisTet == -1)
                                {
                                    FatalErrorInFunction
                                        << "For cell:" << celli
                                        << " at:" << mesh_.cellCentres()[celli]
                                        << " patch:" << patchi
                                        << " face:" << facei
                                        << " at:" << mesh_.faceCentres()[facei]
                                        << " did not find edge:" << e
                                        << " at:" << e.line(mesh_.points())
                                        << " in OWNER-side decomposed faces:"
                                        << flatOutput(subFaces)
                                        << exit(FatalError);
                                }
                            }


                            // Get added cell on otherFacei
                            if (otherCells.size() == 1)
                            {
                                otherTet = otherCells[0];
                            }
                            else
                            {
                                splitFace
                                (
                                    boundaryQuads,
                                    boundaryTris,
                                    otherFacei,
                                    mesh_.boundaryMesh().whichPatch(otherFacei),
                                    quadFaces,  // work space
                                    triFaces,   // work space
                                    subFaces
                                );

                                // Find the edge. Should be only on one of the
                                // subfaces.
                                forAll(subFaces, subFacei)
                                {
                                    const auto& subF = subFaces[subFacei];
                                    if (subF.find(e) != -1)
                                    {
                                        otherTet = otherCells[subFacei];
                                        break;
                                    }
                                }
                                if (otherTet == -1)
                                {
                                    FatalErrorInFunction
                                        << "For cell:" << celli
                                        << " at:" << mesh_.cellCentres()[celli]
                                        << " patch:"
                                        << mesh_.boundaryMesh().whichPatch(otherFacei)
                                        << " face:" << otherFacei
                                        << " at:"
                                        << mesh_.faceCentres()[otherFacei]
                                        << " did not find edge:" << e
                                        << " at:" << e.line(mesh_.points())
                                        << " in decomposed faces:"
                                        << flatOutput(subFaces)
                                        << exit(FatalError);
                                }
                            }
                        }
                        else
                        {
                            thisTet = thisCells[0];
                            otherTet = otherCells[0];
                        }

                        addFace
                        (
                            meshMod,
                            triangle,
                            facei,
                            otherTet,
                            thisTet,
                            -1,         //masterPoint
                            -1,         //fEdges[fp], //masterEdge
                            facei,      //masterFace
                            -1,         //patchi
                            -1,         //zone
                            false
                        );
                    }
                }
            }
        }
    }
}


void Foam::tetDecomposer::updateMesh(const mapPolyMesh& map)
{
    inplaceRenumber(map.reversePointMap(), cellToPoint_);
    inplaceRenumber(map.reversePointMap(), faceToPoint_);

    forAll(faceOwnerCells_, facei)
    {
        inplaceRenumber(map.reverseCellMap(), faceOwnerCells_[facei]);
    }
    forAll(faceNeighbourCells_, facei)
    {
        inplaceRenumber(map.reverseCellMap(), faceNeighbourCells_[facei]);
    }
}


// ************************************************************************* //
