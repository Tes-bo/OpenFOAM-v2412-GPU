/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "pyrMatcher.H"
#include "cellMatcher.H"
#include "primitiveMesh.H"
#include "cellModel.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check (4 tri, 1 quad)
static inline bool checkFaceSizeMatch(const UList<face>& faces)
{
    if (faces.size() != 5)  // facePerCell
    {
        return false;
    }

    int nTris = 0;
    int nQuads = 0;

    for (const face& f : faces)
    {
        const label size = f.size();

        if (size == 3)
        {
            ++nTris;
        }
        else if (size == 4)
        {
            ++nQuads;
        }
        else
        {
            return false;
        }
    }

    return (nTris == 4 && nQuads == 1);
}


// Check (4 tri, 1 quad)
static inline bool checkFaceSizeMatch
(
    const UList<face>& meshFaces,
    const labelUList& cellFaces
)
{
    if (cellFaces.size() != 5)  // facePerCell
    {
        return false;
    }

    int nTris = 0;
    int nQuads = 0;

    for (const label facei : cellFaces)
    {
        const label size = meshFaces[facei].size();

        if (size == 3)
        {
            ++nTris;
        }
        else if (size == 4)
        {
            ++nQuads;
        }
        else
        {
            return false;
        }
    }

    return (nTris == 4 && nQuads == 1);
}

} // End namespace Foam


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

bool Foam::pyrMatcher::test(const UList<face>& faces)
{
    return checkFaceSizeMatch(faces);
}

bool Foam::pyrMatcher::test(const primitiveMesh& mesh, const label celli)
{
    return checkFaceSizeMatch(mesh.faces(), mesh.cells()[celli]);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pyrMatcher::pyrMatcher()
:
    cellMatcher
    (
        vertPerCell,
        facePerCell,
        maxVertPerFace,
        "pyr" // == cellModel::modelNames[cellModel::PYR]
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pyrMatcher::matchShape
(
    const bool checkOnly,
    const faceList& faces,
    const labelList& owner,
    const label celli,
    const labelList& myFaces
)
{
    if (!faceSizeMatch(faces, myFaces))
    {
        return false;
    }

    // Is pyr for sure since no other shape with 1 quad, 4 triangles
    if (checkOnly)
    {
        return true;
    }

    // Calculate localFaces_ and mapping pointMap_, faceMap_
    label numVert = calcLocalFaces(faces, myFaces);

    if (numVert != vertPerCell)
    {
        return false;
    }

    // Set up 'edge' to face mapping.
    calcEdgeAddressing(numVert);

    // Set up point on face to index-in-face mapping
    calcPointFaceIndex();

    // Storage for maps -vertex to mesh and -face to mesh
    vertLabels_.setSize(vertPerCell);
    faceLabels_.setSize(facePerCell);

    //
    // Start from quad face (face0)
    //

    label face0I = -1;
    forAll(faceSize_, facei)
    {
        if (faceSize_[facei] == 4)
        {
            face0I = facei;
            break;
        }
    }
    const face& face0 = localFaces_[face0I];
    label face0vert0 = 0;


    //
    // Try to follow prespecified path on faces of cell,
    // starting at face0vert0
    //

    vertLabels_[0] = pointMap_[face0[face0vert0]];
    faceLabels_[0] = faceMap_[face0I];

    // Walk face 0 from vertex 0 to 1
    label face0vert1 =
        nextVert
        (
            face0vert0,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[1] = pointMap_[face0[face0vert1]];

    // Walk face 0 from vertex 1 to 2
    label face0vert2 =
        nextVert
        (
            face0vert1,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[2] = pointMap_[face0[face0vert2]];

    // Walk face 0 from vertex 2 to 3
    label face0vert3 =
        nextVert
        (
            face0vert2,
            faceSize_[face0I],
            !(owner[faceMap_[face0I]] == celli)
        );
    vertLabels_[3] = pointMap_[face0[face0vert3]];

    // Jump edge from face0 to face1
    label face1I =
        otherFace
        (
            numVert,
            face0[face0vert3],
            face0[face0vert0],
            face0I
        );
    faceLabels_[1] = faceMap_[face1I];

    // Jump edge from face0 to face2
    label face2I =
        otherFace
        (
            numVert,
            face0[face0vert2],
            face0[face0vert3],
            face0I
        );
    faceLabels_[2] = faceMap_[face2I];

    // Jump edge from face0 to face3
    label face3I =
        otherFace
        (
            numVert,
            face0[face0vert1],
            face0[face0vert2],
            face0I
        );
    faceLabels_[3] = faceMap_[face3I];

    // Jump edge from face0 to face4
    label face4I =
        otherFace
        (
            numVert,
            face0[face0vert0],
            face0[face0vert1],
            face0I
        );
    faceLabels_[4] = faceMap_[face4I];

    const face& face4 = localFaces_[face4I];

    // Get index of vert0 in face 4
    label face4vert0 = pointFaceIndex_[face0[face0vert0]][face4I];

    // Walk face 4 from vertex 0 to 4
    label face4vert4 =
        nextVert
        (
            face4vert0,
            faceSize_[face4I],
            !(owner[faceMap_[face4I]] == celli)
        );
    vertLabels_[4] = pointMap_[face4[face4vert4]];

    return true;
}


Foam::label Foam::pyrMatcher::faceHashValue() const
{
    return 4*3+4;
}


bool Foam::pyrMatcher::faceSizeMatch
(
    const faceList& meshFaces,
    const labelList& cellFaces
) const
{
    return checkFaceSizeMatch(meshFaces, cellFaces);
}


bool Foam::pyrMatcher::matches
(
    const primitiveMesh& mesh,
    const label celli,
    cellShape& shape
)
{
    if
    (
        matchShape
        (
            false,
            mesh.faces(),
            mesh.faceOwner(),
            celli,
            mesh.cells()[celli]
        )
    )
    {
        shape.reset(model(), vertLabels());
        return true;
    }

    return false;
}


// ************************************************************************* //
