/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

Description
    Create the list of loops of outside vertices. Goes wrong on multiply
    connected edges (loops will be unclosed).

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcEdgeLoops() const
{
    DebugInFunction << "Calculating boundary edge loops" << endl;

    if (edgeLoopsPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "edge loops already calculated"
            << abort(FatalError);
    }

    const edgeList& patchEdges = edges();
    const label nIntEdges = nInternalEdges();
    const label nBdryEdges = (patchEdges.size() - nIntEdges);

    // Size return list plenty big
    edgeLoopsPtr_.reset(new labelListList(nBdryEdges));
    auto& edgeLoops = *edgeLoopsPtr_;

    if (nBdryEdges == 0)
    {
        return;
    }

    const labelListList& patchPointEdges = pointEdges();


    //
    // Walk point-edge-point and assign loop number
    //

    // Temporary storage for vertices of current loop
    DynamicList<label> loop(nBdryEdges);

    // In a loop? - per boundary edge
    boolList unvisited(nBdryEdges, true);

    // Current loop number
    label numLoops = 0;

    // Walk all boundary edges not yet in a loop
    for
    (
        label bndEdgei = -1;
        (bndEdgei = unvisited.find(true)) >= 0;
        /*nil*/
    )
    {
        label currentEdgei = (bndEdgei + nIntEdges);

        // Walk from first all the way round, assigning loops
        label currentVerti = patchEdges[currentEdgei].first();

        loop.clear();

        do
        {
            loop.push_back(currentVerti);

            unvisited[currentEdgei - nIntEdges] = false;

            // Step to next vertex
            currentVerti = patchEdges[currentEdgei].otherVertex(currentVerti);

            // Step to next (unmarked, boundary) edge.
            currentEdgei = -1;

            for (const label edgei : patchPointEdges[currentVerti])
            {
                if (edgei >= nIntEdges && unvisited[edgei - nIntEdges])
                {
                    // Unvisited boundary edge
                    currentEdgei = edgei;
                    break;
                }
            }
        }
        while (currentEdgei != -1);

        // Done all for current loop - copy to edgeLoops
        edgeLoops[numLoops] = loop;

        ++numLoops;
    }

    edgeLoops.resize(numLoops);

    DebugInFunction << "Calculated boundary edge loops" << nl;
}


template<class FaceList, class PointField>
const Foam::labelListList&
Foam::PrimitivePatch<FaceList, PointField>::edgeLoops() const
{
    if (!edgeLoopsPtr_)
    {
        calcEdgeLoops();
    }

    return *edgeLoopsPtr_;
}


// ************************************************************************* //
