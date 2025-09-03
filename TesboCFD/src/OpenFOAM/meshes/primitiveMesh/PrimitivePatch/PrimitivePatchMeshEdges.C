/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
Foam::edge
Foam::PrimitivePatch<FaceList, PointField>::meshEdge(const label edgei) const
{
    return Foam::edge(this->meshPoints(), this->edges()[edgei]);
}


template<class FaceList, class PointField>
Foam::edge
Foam::PrimitivePatch<FaceList, PointField>::meshEdge(const edge& e) const
{
    return Foam::edge(this->meshPoints(), e);
}


template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::
meshEdges
(
    const edgeList& allEdges,
    const labelListList& cellEdges,
    const labelList& faceCells
) const
{
    DebugInFunction
        << "Calculating labels of patch edges in mesh edge list" << nl;

    // The output storage
    labelList meshEdgeLabels(this->nEdges());

    const labelListList& EdgeFaces = edgeFaces();

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll(meshEdgeLabels, edgei)
    {
        const edge globalEdge(this->meshEdge(edgei));

        bool found = false;

        // For each patch face sharing the edge
        for (const label patchFacei : EdgeFaces[edgei])
        {
            // The cell next to the face
            const label curCelli = faceCells[patchFacei];

            // Check the cell edges
            for (const label cellEdgei : cellEdges[curCelli])
            {
                if (allEdges[cellEdgei] == globalEdge)
                {
                    found = true;
                    meshEdgeLabels[edgei] = cellEdgei;
                    break;
                }
            }

            if (found) break;
        }
    }

    return meshEdgeLabels;
}


template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::meshEdges
(
    const edgeList& allEdges,
    const labelListList& pointEdges
) const
{
    DebugInFunction
        << "Calculating labels of patch edges in mesh edge list" << nl;

    labelList meshEdgeLabels(this->nEdges());

    // WARNING: Remember that local edges address into local point list;
    // local-to-global point label translation is necessary
    forAll(meshEdgeLabels, edgei)
    {
        const edge globalEdge(this->meshEdge(edgei));

        // Check the attached edges
        for (const label meshEdgei : pointEdges[globalEdge.start()])
        {
            if (allEdges[meshEdgei] == globalEdge)
            {
                meshEdgeLabels[edgei] = meshEdgei;
                break;
            }
        }
    }

    return meshEdgeLabels;
}


template<class FaceList, class PointField>
Foam::label
Foam::PrimitivePatch<FaceList, PointField>::meshEdge
(
    const label edgei,
    const edgeList& allEdges,
    const labelListList& pointEdges
) const
{
    // Need local-to-global point label translation
    const edge globalEdge(this->meshEdge(edgei));

    // Check attached edges
    for (const label meshEdgei : pointEdges[globalEdge.start()])
    {
        if (allEdges[meshEdgei] == globalEdge)
        {
            return meshEdgei;
        }
    }

    return -1;
}


template<class FaceList, class PointField>
Foam::labelList
Foam::PrimitivePatch<FaceList, PointField>::meshEdges
(
    const labelUList& edgeLabels,
    const edgeList& allEdges,
    const labelListList& pointEdges
) const
{
    labelList meshEdgeLabels(edgeLabels.size());

    forAll(meshEdgeLabels, edgei)
    {
        meshEdgeLabels[edgei] =
            this->meshEdge(edgeLabels[edgei], allEdges, pointEdges);
    }

    return meshEdgeLabels;
}


template<class FaceList, class PointField>
Foam::label
Foam::PrimitivePatch<FaceList, PointField>::findEdge
(
    const edge& e
) const
{
    if (e.good() && e.first() < nPoints() && e.second() < nPoints())
    {
        // Get pointEdges from the starting point and search all the candidates
        const edgeList& myEdges = this->edges();

        for (const label patchEdgei : pointEdges()[e.first()])
        {
            if (e == myEdges[patchEdgei])
            {
                return patchEdgei;
            }
        }
    }

    return -1;  // Not found, or invalid edge
}


// ************************************************************************* //
