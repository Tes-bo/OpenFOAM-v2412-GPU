/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2018 OpenCFD Ltd.
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
    Searching and marking zones of the patch.

\*---------------------------------------------------------------------------*/

#include "PatchTools.H"
#include "bitSet.H"
#include "boundBox.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class BoolListType, class FaceList, class PointField>
void Foam::PatchTools::markZone
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    const label facei,
    const label currentZone,
    labelList&  faceZone
)
{
    const labelListList& faceEdges = p.faceEdges();
    const labelListList& edgeFaces = p.edgeFaces();

    // List of faces whose faceZone has been set.
    labelList changedFaces(1, facei);

    while (true)
    {
        // Pick up neighbours of changedFaces
        DynamicList<label> newChangedFaces(2*changedFaces.size());

        forAll(changedFaces, i)
        {
            label facei = changedFaces[i];

            const labelList& fEdges = faceEdges[facei];

            forAll(fEdges, fEdgeI)
            {
                label edgeI = fEdges[fEdgeI];

                if (!borderEdge[edgeI])
                {
                    const labelList& eFaceLst = edgeFaces[edgeI];

                    forAll(eFaceLst, j)
                    {
                        label nbrFacei = eFaceLst[j];

                        if (faceZone[nbrFacei] == -1)
                        {
                            faceZone[nbrFacei] = currentZone;
                            newChangedFaces.append(nbrFacei);
                        }
                        else if (faceZone[nbrFacei] != currentZone)
                        {
                            FatalErrorInFunction
                                << "Zones " << faceZone[nbrFacei]
                                << " at face " << nbrFacei
                                << " connects to zone " << currentZone
                                << " at face " << facei
                                << abort(FatalError);
                        }
                    }
                }
            }
        }

        if (newChangedFaces.empty())
        {
            break;
        }

        // transfer from dynamic to normal list
        changedFaces.transfer(newChangedFaces);
    }
}


template<class BoolListType, class FaceList, class PointField>
Foam::label
Foam::PatchTools::markZones
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& borderEdge,
    labelList& faceZone
)
{
    faceZone.setSize(p.size());
    faceZone = -1;

    label zoneI = 0;
    for (label startFacei = 0; startFacei < faceZone.size();)
    {
        // Find next non-visited face
        for (; startFacei < faceZone.size(); ++startFacei)
        {
            if (faceZone[startFacei] == -1)
            {
                faceZone[startFacei] = zoneI;
                markZone(p, borderEdge, startFacei, zoneI, faceZone);
                zoneI++;
                break;
            }
        }
    }

    return zoneI;
}


template<class BoolListType, class FaceList, class PointField>
void
Foam::PatchTools::subsetMap
(
    const PrimitivePatch<FaceList, PointField>& p,
    const BoolListType& includeFaces,
    labelList& pointMap,
    labelList& faceMap
)
{
    const auto& localFaces = p.localFaces();

    faceMap.resize(localFaces.size());
    pointMap.clear();

    bitSet pointUsed(p.nPoints());

    label facei  = 0;

    forAll(localFaces, oldFacei)
    {
        if (includeFaces[oldFacei])
        {
            // Compact storage for new faces
            faceMap[facei++] = oldFacei;

            // Local points used by face
            pointUsed.set(localFaces[oldFacei]);
        }
    }

    // The newToOld mappings
    faceMap.resize(facei);
    pointMap = pointUsed.sortedToc();
}


template<class FaceList, class PointField>
void Foam::PatchTools::calcBounds
(
    const PrimitivePatch<FaceList, PointField>& p,
    boundBox& bb,
    label& nPoints
)
{
    // Unfortunately nPoints constructs meshPoints() so do compact version
    // ourselves
    const PointField& points = p.points();

    bitSet pointUsed(points.size());

    nPoints = 0;
    bb.reset();

    for (const auto& f : p)
    {
        for (const label pointi : f)
        {
            if (pointUsed.set(pointi))
            {
                bb.add(points[pointi]);
                ++nPoints;
            }
        }
    }
}


// ************************************************************************* //
