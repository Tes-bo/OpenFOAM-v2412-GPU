/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "patchEdgeSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "searchableSurface.H"
#include "Time.H"
#include "mergePoints.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchEdgeSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchEdgeSet, word);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchEdgeSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    const labelList patchIDs(patchSet_.sortedToc());

    for (const label patchi : patchIDs)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        const edgeList& edges = pp.edges();
        const labelList& mp = pp.meshPoints();
        const pointField& pts = pp.points();

        pointField start(edges.size());
        pointField end(edges.size());
        forAll(edges, edgei)
        {
            const edge& e = edges[edgei];
            start[edgei] = pts[mp[e[0]]];
            end[edgei] = pts[mp[e[1]]];
        }


        List<List<pointIndexHit>> hits;
        surfPtr_->findLineAll(start, end, hits);

        forAll(hits, edgei)
        {
            const List<pointIndexHit>& eHits = hits[edgei];

            if (eHits.size())
            {
                const label meshFacei = pp.start()+pp.edgeFaces()[edgei][0];
                const label celli = mesh().faceOwner()[meshFacei];
                for (const auto& hit : eHits)
                {
                    const point& pt = hit.hitPoint();
                    samplingPts.append(pt);
                    samplingCells.append(celli);
                    samplingFaces.append(meshFacei);
                    samplingSegments.append(0);
                    samplingCurveDist.append(mag(pt-origin_));
                }
            }
        }

        ////- Alternative hardcoded to use plane
        //forAll(edges, edgei)
        //{
        //    const edge& e = edges[edgei];
        //    const point& p0 = pts[mp[e[0]]];
        //    const point& p1 = pts[mp[e[1]]];
        //    const vector dir(p1-p0);
        //    const scalar s = pl_.normalIntersect(p0, dir);
        //
        //    if (s >= 0.0 && s < 1.0)
        //    {
        //        const point pt(p0+s*dir);
        //
        //        // Calculate distance on curve
        //        //const scalar dist = (pt-start_)&curve;
        //        //if (dist >= 0 && dist < magCurve)
        //        const scalar dist = mag(pt-pl_.origin());
        //
        //        // Take any face using this edge
        //        const label meshFacei = pp.start()+pp.edgeFaces()[edgei][0];
        //
        //        samplingPts.append(pt);
        //        samplingCells.append(mesh().faceOwner()[meshFacei]);
        //        samplingFaces.append(meshFacei);
        //        samplingSegments.append(0);
        //        samplingCurveDist.append(dist);
        //    }
        //}
    }

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();


    labelList pointMap;
    const label nMerged = Foam::mergePoints
    (
        samplingPts,
        SMALL,          // mergeTol
        false,          // verbose = false
        pointMap
    );

    if (nMerged == samplingPts.size())
    {
        // Nothing merged
        setSamples
        (
            samplingPts,
            samplingCells,
            samplingFaces,
            samplingSegments,
            samplingCurveDist
        );
    }
    else
    {
        // Compress out duplicates

        List<point> newSamplingPts(nMerged);
        List<label> newSamplingCells(nMerged);
        List<label> newSamplingFaces(nMerged);
        List<label> newSamplingSegments(nMerged);
        List<scalar> newSamplingCurveDist(nMerged);

        forAll(pointMap, i)
        {
            const label newi = pointMap[i];
            newSamplingPts[newi] = samplingPts[i];
            newSamplingCells[newi] = samplingCells[i];
            newSamplingFaces[newi] = samplingFaces[i];
            newSamplingSegments[newi] = samplingSegments[i];
            newSamplingCurveDist[newi] = samplingCurveDist[i];
        }

        setSamples
        (
            newSamplingPts,
            newSamplingCells,
            newSamplingFaces,
            newSamplingSegments,
            newSamplingCurveDist
        );
    }

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchEdgeSet::patchEdgeSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    surfPtr_
    (
        searchableSurface::New
        (
            dict.get<word>("surfaceType"),
            IOobject
            (
                dict.getOrDefault("surfaceName", name),
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    origin_(dict.get<point>("origin")),
    //end_(dict.get<point>("end")),
    //pl_(dict),
    patchSet_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    )
{
    genSamples();
}


// ************************************************************************* //
