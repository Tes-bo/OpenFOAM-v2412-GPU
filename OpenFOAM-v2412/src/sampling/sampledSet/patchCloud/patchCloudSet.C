/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "patchCloudSet.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "Time.H"
#include "meshTools.H"
// For 'nearInfo' helper class only
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchCloudSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchCloudSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchCloudSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
) const
{
    DebugInfo << "patchCloudSet : sampling on patches :" << endl;

    // Construct search tree for all patch faces.
    label sz = 0;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        sz += pp.size();

        DebugInfo << "    " << pp.name() << " size " << pp.size() << endl;
    }

    labelList patchFaces(sz);
    sz = 0;
    treeBoundBox bb;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];

        forAll(pp, i)
        {
            patchFaces[sz++] = pp.start()+i;
        }

        // Without reduction.
        bb.add(pp.points(), pp.meshPoints());
    }

    // Not very random
    Random rndGen(123456);
    // Make bb asymmetric just to avoid problems on symmetric meshes
    // Make sure bb is 3D.
    bb.inflate(rndGen, 1e-4, ROOTVSMALL);


    indexedOctree<treeDataFace> patchTree
    (
        treeDataFace(mesh(), patchFaces),  // boundary faces only

        bb,             // overall search domain
        8,              // maxLevel
        10,             // leafsize
        3.0             // duplicity
    );

    // Force calculation of face-diagonal decomposition
    (void)mesh().tetBasePtIs();


    // All the info for nearest. Construct to miss
    List<mappedPatchBase::nearInfo> nearest(sampleCoords_.size());

    forAll(sampleCoords_, sampleI)
    {
        const auto& treeData = patchTree.shapes();
        const point& sample = sampleCoords_[sampleI];

        pointIndexHit& nearInfo = nearest[sampleI].first();
        auto& distSqrProc = nearest[sampleI].second();

        // Find the nearest locally
        if (patchFaces.size())
        {
            nearInfo = patchTree.findNearest(sample, sqr(searchDist_));
        }
        else
        {
            nearInfo.setMiss();
        }


        // Fill in the distance field and the processor field
        if (!nearInfo.hit())
        {
            distSqrProc.first() = Foam::sqr(GREAT);
            distSqrProc.second() = Pstream::myProcNo();
        }
        else
        {
            // Set nearest to mesh face label
            const label objectIndex = treeData.objectIndex(nearInfo.index());

            nearInfo.setIndex(objectIndex);

            distSqrProc.first() = sample.distSqr(nearInfo.point());
            distSqrProc.second() = Pstream::myProcNo();
        }
    }


    // Find nearest - globally consistent
    Pstream::listCombineReduce(nearest, mappedPatchBase::nearestEqOp());


    if (debug && Pstream::master())
    {
        OFstream str
        (
            mesh().time().path()
          / name()
          + "_nearest.obj"
        );
        Info<< "Dumping mapping as lines from supplied points to"
            << " nearest patch face to file " << str.name() << endl;

        label vertI = 0;

        forAll(nearest, i)
        {
            if (nearest[i].first().hit())
            {
                meshTools::writeOBJ(str, sampleCoords_[i]);
                ++vertI;
                meshTools::writeOBJ(str, nearest[i].first().point());
                ++vertI;
                str << "l " << vertI-1 << ' ' << vertI << nl;
            }
        }
    }


    // Store the sampling locations on the nearest processor
    forAll(nearest, sampleI)
    {
        const pointIndexHit& nearInfo = nearest[sampleI].first();

        if (nearInfo.hit())
        {
            if (nearest[sampleI].second().second() == Pstream::myProcNo())
            {
                label facei = nearInfo.index();

                samplingPts.append(nearInfo.point());
                samplingCells.append(mesh().faceOwner()[facei]);
                samplingFaces.append(facei);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0 * sampleI);
            }
        }
        else
        {
            // No processor found point near enough. Mark with special value
            // which is intercepted when interpolating
            if (Pstream::master())
            {
                samplingPts.append(sampleCoords_[sampleI]);
                samplingCells.append(-1);
                samplingFaces.append(-1);
                samplingSegments.append(0);
                samplingCurveDist.append(1.0 * sampleI);
            }
        }
    }
}


void Foam::patchCloudSet::genSamples()
{
    // Storage for sample points
    DynamicList<point> samplingPts;
    DynamicList<label> samplingCells;
    DynamicList<label> samplingFaces;
    DynamicList<label> samplingSegments;
    DynamicList<scalar> samplingCurveDist;

    calcSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    samplingPts.shrink();
    samplingCells.shrink();
    samplingFaces.shrink();
    samplingSegments.shrink();
    samplingCurveDist.shrink();

    setSamples
    (
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    if (debug)
    {
        write(Info);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const word& axis,
    const List<point>& sampleCoords,
    const labelHashSet& patchSet,
    const scalar searchDist
)
:
    sampledSet(name, mesh, searchEngine, axis),
    sampleCoords_(sampleCoords),
    patchSet_(patchSet),
    searchDist_(searchDist)
{
    genSamples();
}


Foam::patchCloudSet::patchCloudSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    sampleCoords_(dict.get<pointField>("points")),
    patchSet_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    ),
    searchDist_(dict.get<scalar>("maxDistance"))
{
    genSamples();
}


// ************************************************************************* //
