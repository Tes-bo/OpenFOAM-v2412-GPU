/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "PatchTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template
<
    class FaceList1, class PointField1,
    class FaceList2, class PointField2
>
void Foam::PatchTools::matchPoints
(
    const PrimitivePatch<FaceList1, PointField1>& p1,
    const PrimitivePatch<FaceList2, PointField2>& p2,

    labelList& p1PointLabels,
    labelList& p2PointLabels
)
{
    p1PointLabels.resize(p1.nPoints());
    p2PointLabels.resize(p1.nPoints());

    label nMatches = 0;

    forAll(p1.meshPoints(), pointi)
    {
        const label meshPointi = p1.meshPoints()[pointi];

        const auto iter = p2.meshPointMap().cfind(meshPointi);

        if (iter.good())
        {
            p1PointLabels[nMatches] = pointi;
            p2PointLabels[nMatches] = iter.val();
            ++nMatches;
        }
    }

    p1PointLabels.resize(nMatches);
    p2PointLabels.resize(nMatches);
}


template
<
    class FaceList1, class PointField1,
    class FaceList2, class PointField2
>
void Foam::PatchTools::matchEdges
(
    const PrimitivePatch<FaceList1, PointField1>& p1,
    const PrimitivePatch<FaceList2, PointField2>& p2,

    labelList& p1EdgeLabels,
    labelList& p2EdgeLabels,
    bitSet& sameOrientation
)
{
    p1EdgeLabels.resize(p1.nEdges());
    p2EdgeLabels.resize(p1.nEdges());
    sameOrientation.resize(p1.nEdges());
    sameOrientation = false;

    label nMatches = 0;

    EdgeMap<label> edgeToIndex(2*p1.nEdges());
    forAll(p1.edges(), edgei)
    {
        // Map lookup with globalEdge
        edgeToIndex.insert(p1.meshEdge(edgei), edgei);
    }

    forAll(p2.edges(), edgei)
    {
        const edge meshEdge2(p2.meshEdge(edgei));

        const auto iter = edgeToIndex.cfind(meshEdge2);

        if (iter.good())
        {
            p1EdgeLabels[nMatches] = iter.val();
            p2EdgeLabels[nMatches] = edgei;
            sameOrientation.set(nMatches, (meshEdge2[0] == iter.key()[0]));
            ++nMatches;
        }
    }

    p1EdgeLabels.resize(nMatches);
    p2EdgeLabels.resize(nMatches);
    sameOrientation.resize(nMatches);
}


// ************************************************************************* //
