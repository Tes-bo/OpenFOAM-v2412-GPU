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
    Checks topology of the patch.

\*---------------------------------------------------------------------------*/

#include "PrimitivePatch.H"
#include "Map.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class FaceList, class PointField>
void
Foam::PrimitivePatch<FaceList, PointField>::visitPointRegion
(
    const label pointi,
    const labelUList& pFaces,
    const label startFacei,
    const label startEdgei,
    UList<bool>& pFacesVisited
) const
{
    const label index = pFaces.find(startFacei);

    if (index >= 0 && !pFacesVisited[index])
    {
        // Mark face as visited
        pFacesVisited[index] = true;

        // Step to next edge on face which is still using pointi
        label nextEdgei = -1;

        for (const label faceEdgei : faceEdges()[startFacei])
        {
            const edge& e = edges()[faceEdgei];

            if (faceEdgei != startEdgei && e.contains(pointi))
            {
                nextEdgei = faceEdgei;
                break;
            }
        }

        if (nextEdgei == -1)
        {
            FatalErrorInFunction
                << "Problem: cannot find edge out of "
                << faceEdges()[startFacei]
                << "on face " << startFacei << " that uses point " << pointi
                << " and is not edge " << startEdgei << abort(FatalError);
        }

        // Walk to next face(s) across edge.
        for (const label edgeFacei : edgeFaces()[nextEdgei])
        {
            if (edgeFacei != startFacei)
            {
                visitPointRegion
                (
                    pointi,
                    pFaces,
                    edgeFacei,
                    nextEdgei,
                    pFacesVisited
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class FaceList, class PointField>
typename Foam::PrimitivePatch<FaceList, PointField>::surfaceTopo
Foam::PrimitivePatch<FaceList, PointField>::surfaceType
(
    labelHashSet* badEdgesPtr
) const
{
    if (badEdgesPtr)
    {
        badEdgesPtr->clear();
    }

    bool foundError = false;

    surfaceTopo pType = surfaceTopo::MANIFOLD;

    const labelListList& edgeFcs = edgeFaces();

    forAll(edgeFcs, edgei)
    {
        const label nNbrs = edgeFcs[edgei].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            // Surface is illegal
            foundError = true;

            if (badEdgesPtr)
            {
                // Record and continue
                if (nNbrs > 2)
                {
                    badEdgesPtr->insert(edgei);
                }
            }
            else
            {
                // Early exit when not recording bad edges
                break;
            }
        }
        else if (nNbrs == 1)
        {
            // Surface might be open or illegal so keep looping.
            pType = surfaceTopo::OPEN;
        }
    }

    if (foundError)
    {
        return surfaceTopo::ILLEGAL;
    }

    return pType;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::checkTopology
(
    const bool report,
    labelHashSet* pointSetPtr
) const
{
    // Check edgeFaces

    bool foundError = false;

    const labelListList& edgeFcs = edgeFaces();

    forAll(edgeFcs, edgei)
    {
        const label nNbrs = edgeFcs[edgei].size();

        if (nNbrs < 1 || nNbrs > 2)
        {
            foundError = true;

            if (report)
            {
                Info<< "Edge " << edgei << " with vertices:" << edges()[edgei]
                    << " has " << nNbrs << " face neighbours" << endl;
            }

            if (pointSetPtr)
            {
                const edge& e = edges()[edgei];

                pointSetPtr->insert(meshPoints()[e.first()]);
                pointSetPtr->insert(meshPoints()[e.second()]);
            }
        }
    }

    return foundError;
}


template<class FaceList, class PointField>
bool
Foam::PrimitivePatch<FaceList, PointField>::checkPointManifold
(
    const bool report,
    labelHashSet* pointSetPtr
) const
{
    const labelListList& pf = pointFaces();
    const labelListList& pe = pointEdges();
    const labelListList& ef = edgeFaces();
    const labelList& mp = meshPoints();

    bool foundError = false;

    // Visited faces (as indices into pFaces)
    DynamicList<bool> pFacesVisited;

    forAll(pf, pointi)
    {
        const labelList& pFaces = pf[pointi];

        pFacesVisited.resize_nocopy(pFaces.size());
        pFacesVisited = false;

        // Starting edge
        const labelList& pEdges = pe[pointi];
        const label startEdgei = pEdges[0];

        const labelList& eFaces = ef[startEdgei];

        for (const label edgeFacei : eFaces)
        {
            // Visit all faces using pointi, starting from eFaces[i] and
            // startEdgei. Mark off all faces visited in pFacesVisited.
            this->visitPointRegion
            (
                pointi,
                pFaces,
                edgeFacei,  // starting face for walk
                startEdgei, // starting edge for walk
                pFacesVisited
            );
        }

        // After this all faces using pointi should have been visited and
        // marked off in pFacesVisited.

        if (pFacesVisited.contains(false))
        {
            foundError = true;

            const label meshPointi = mp[pointi];

            if (pointSetPtr)
            {
                pointSetPtr->insert(meshPointi);
            }

            if (report)
            {
                Info<< "Point " << meshPointi
                    << " uses faces which are not connected through an edge"
                    << nl
                    << "This means that the surface formed by this patched"
                    << " is multiply connected at this point" << nl
                    << "Connected (patch) faces:" << nl;

                forAll(pFacesVisited, i)
                {
                    if (pFacesVisited[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }

                Info<< nl << "Unconnected (patch) faces:" << nl;
                forAll(pFacesVisited, i)
                {
                    if (!pFacesVisited[i])
                    {
                        Info<< "    " << pFaces[i] << endl;
                    }
                }
            }
        }
    }

    return foundError;
}


// ************************************************************************* //
