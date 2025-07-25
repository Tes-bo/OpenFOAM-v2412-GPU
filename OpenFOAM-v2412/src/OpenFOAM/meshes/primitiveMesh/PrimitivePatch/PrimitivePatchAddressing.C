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

Description
    This function calculates the list of patch edges, defined on the list of
    points supporting the patch. The edges are ordered:
    - 0..nInternalEdges-1 : upper triangular order
    - nInternalEdges..    : boundary edges (no particular order)

    Other patch addressing information is also calculated:
    - faceFaces with neighbour faces in ascending order
    - edgeFaces with ascending face order
    - faceEdges sorted according to edges of a face

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::calcAddressing() const
{
    DebugInFunction << "Calculating patch addressing" << nl;

    if (edgesPtr_ || faceFacesPtr_ || edgeFacesPtr_ || faceEdgesPtr_)
    {
        // An error to recalculate if already allocated
        FatalErrorInFunction
            << "addressing already calculated"
            << abort(FatalError);
    }

    // get reference to localFaces
    const List<face_type>& locFcs = localFaces();

    // get reference to pointFaces
    const labelListList& pf = pointFaces();

    // Guess the max number of edges/neighbours for a face
    label maxEdges = 0;
    for (const auto& f : locFcs)
    {
        maxEdges += f.nEdges();
    }

    // create the lists for the various results. (resized on completion)
    edgesPtr_.reset(new edgeList(maxEdges));
    auto& edges = *edgesPtr_;

    edgeFacesPtr_.reset(new labelListList(maxEdges));
    auto& edgeFaces = *edgeFacesPtr_;

    faceEdgesPtr_.reset(new labelListList(locFcs.size()));
    auto& faceEdges = *faceEdgesPtr_;

    // initialise the lists of subshapes for each face to avoid duplication
    edgeListList faceIntoEdges(locFcs.size());

    forAll(locFcs, facei)
    {
        faceIntoEdges[facei] = locFcs[facei].edges();
        faceEdges[facei].resize(faceIntoEdges[facei].size(), -1);
    }

    // faceFaces created using a dynamic list.  Cannot guess size because
    // of multiple connections
    List<DynamicList<label>> ff(locFcs.size());


    // This algorithm will produce a separated list of edges, internal edges
    // starting from 0 and boundary edges starting from the top and
    // growing down.

    label nEdges = 0;

    // Note that faceIntoEdges is sorted acc. to local vertex numbering
    // in face (i.e. curEdges[0] is edge between f[0] and f[1])

    // For all local faces ...
    forAll(locFcs, facei)
    {
        // Get reference to vertices of current face and corresponding edges.
        const face_type& curF = locFcs[facei];
        const edgeList& curEdges = faceIntoEdges[facei];

        // Record the neighbour face.  Multiple connectivity allowed
        List<DynamicList<label>> neiFaces(curF.size());
        List<DynamicList<label>> edgeOfNeiFace(curF.size());

        label nNeighbours = 0;

        // For all edges ...
        forAll(curEdges, edgeI)
        {
            // If the edge is already detected, skip
            if (faceEdges[facei][edgeI] >= 0) continue;

            // Set reference to the current edge
            const edge& e = curEdges[edgeI];

            // Collect neighbours for the current face vertex.

            const labelList& nbrFaces = pf[e.start()];

            bool found = false;

            forAll(nbrFaces, nbrFacei)
            {
                // set reference to the current neighbour
                label curNei = nbrFaces[nbrFacei];

                // Reject neighbours with the lower label
                if (curNei > facei)
                {
                    // get the reference to subshapes of the neighbour
                    const edgeList& searchEdges = faceIntoEdges[curNei];

                    forAll(searchEdges, neiEdgeI)
                    {
                        if (searchEdges[neiEdgeI] == e)
                        {
                            // Match
                            found = true;

                            neiFaces[edgeI].append(curNei);
                            edgeOfNeiFace[edgeI].append(neiEdgeI);

                            // Record faceFaces both ways
                            ff[facei].append(curNei);
                            ff[curNei].append(facei);

                            // Keep searching due to multiple connectivity
                        }
                    }
                }
            } // End of neighbouring faces

            if (found)
            {
                // Register another detected internal edge
                nNeighbours++;
            }
        } // End of current edges

        // Add the edges in increasing number of neighbours.
        // Note: for multiply connected surfaces, the lower index neighbour for
        // an edge will come first.

        // Add the faces in the increasing order of neighbours
        for (label neiSearch = 0; neiSearch < nNeighbours; neiSearch++)
        {
            // Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = locFcs.size();

            forAll(neiFaces, nfI)
            {
                if (neiFaces[nfI].size() && neiFaces[nfI][0] < minNei)
                {
                    nextNei = nfI;
                    minNei = neiFaces[nfI][0];
                }
            }

            if (nextNei > -1)
            {
                // Add the face to the list of faces
                edges[nEdges] = curEdges[nextNei];

                // Set face-edge and face-neighbour-edge to current face label
                faceEdges[facei][nextNei] = nEdges;

                DynamicList<label>& cnf = neiFaces[nextNei];
                DynamicList<label>& eonf = edgeOfNeiFace[nextNei];

                // Set edge-face addressing
                labelList& curEf = edgeFaces[nEdges];
                curEf.resize(cnf.size() + 1);
                curEf[0] = facei;

                forAll(cnf, cnfI)
                {
                    faceEdges[cnf[cnfI]][eonf[cnfI]] = nEdges;

                    curEf[cnfI + 1] = cnf[cnfI];
                }

                // Stop the neighbour from being used again
                cnf.clear();
                eonf.clear();

                // Increment number of faces counter
                nEdges++;
            }
            else
            {
                FatalErrorInFunction
                    << "Error in internal edge insertion"
                    << abort(FatalError);
            }
        }
    }

    nInternalEdges_ = nEdges;

    // Do boundary faces

    forAll(faceEdges, facei)
    {
        labelList& curEdges = faceEdges[facei];

        forAll(curEdges, edgeI)
        {
            if (curEdges[edgeI] < 0)
            {
                // Grab edge and faceEdge
                edges[nEdges] = faceIntoEdges[facei][edgeI];
                curEdges[edgeI] = nEdges;

                // Add edgeFace
                labelList& curEf = edgeFaces[nEdges];
                curEf.resize(1);
                curEf[0] = facei;

                nEdges++;
            }
        }
    }

    // edges
    edges.resize(nEdges);

    // edgeFaces list
    edgeFaces.resize(nEdges);

    // faceFaces list
    faceFacesPtr_.reset(new labelListList(locFcs.size()));
    auto& faceFaces = *faceFacesPtr_;

    forAll(faceFaces, facei)
    {
        faceFaces[facei].transfer(ff[facei]);
    }


    DebugInFunction << "Calculated patch addressing" << nl;
}


// ************************************************************************* //
