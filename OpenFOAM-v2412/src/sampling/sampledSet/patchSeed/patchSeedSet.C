/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "patchSeedSet.H"
#include "polyMesh.H"
#include "treeBoundBox.H"
#include "treeDataFace.H"
#include "mappedPatchBase.H"
#include "indirectPrimitivePatch.H"
#include "triangulatedPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchSeedSet, 0);
    addToRunTimeSelectionTable(sampledSet, patchSeedSet, word);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::patchSeedSet::calcPatchSamples
(
    const label nAvailable,
    const label nPatchPoints,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    if (nAvailable < 1)
    {
        return;
    }

    Random& rndGen = *rndGenPtr_;

    globalIndex globalSampleNumbers(nAvailable);
    label nGlobalPatchPoints = returnReduce(nPatchPoints, sumOp<label>());

    point pt;
    label facei;
    label celli;

    const bool perturb = true;

    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        triangulatedPatch tp(pp, perturb);

        const label np = nAvailable*pp.size()/scalar(nGlobalPatchPoints);
        for (label i = 0; i < np; ++i)
        {
            tp.randomLocalPoint(rndGen, pt, facei, celli);

            samplingPts.append(pt);
            samplingCells.append(celli);
            samplingFaces.append(pp.start() + facei);
            samplingSegments.append(0);
            samplingCurveDist.append(globalSampleNumbers.toGlobal(i));
        }
    }
}


void Foam::patchSeedSet::calcSelectedLocations
(
    const label nAvailable,
    const label nPatchPoints,
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    if (nAvailable < 1)
    {
        return;
    }

    Random& rndGen = *rndGenPtr_;

    labelList patchFaces(nPatchPoints);
    label sz = 0;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        forAll(pp, localFacei)
        {
            patchFaces[sz++] = pp.start() + localFacei;
        }
    }

    {
        DynamicList<label> newPatchFaces(patchFaces.size());

        // Find the nearest patch face
        {
            // 1. All processors find nearest local patch face for all
            //    selectedLocations

            // All the info for nearest. Construct to miss
            List<mappedPatchBase::nearInfo> nearest(nAvailable);

            const indirectPrimitivePatch pp
            (
                IndirectList<face>(mesh().faces(), patchFaces),
                mesh().points()
            );

            treeBoundBox patchBb
            (
                treeBoundBox(pp.points(), pp.meshPoints())
                    .extend(rndGen, 1e-4, ROOTVSMALL)
            );

            indexedOctree<treeDataFace> boundaryTree
            (
                treeDataFace(mesh(), patchFaces),  // boundary faces only

                patchBb,        // overall search domain
                8,              // maxLevel
                10,             // leafsize
                3.0             // duplicity
            );

            // Get some global dimension so all points are equally likely
            // to be found
            const scalar globalDistSqr
            (
                //boundBox(pp.points(), pp.meshPoints(), true).magSqr()
                GREAT
            );

            for (label sampleI = 0; sampleI < nAvailable; ++sampleI)
            {
                const auto& treeData = boundaryTree.shapes();
                const point& sample = selectedLocations_[sampleI];

                pointIndexHit& nearInfo = nearest[sampleI].first();
                auto& distSqrProc = nearest[sampleI].second();

                nearInfo = boundaryTree.findNearest
                (
                    sample,
                    globalDistSqr
                );

                if (!nearInfo.hit())
                {
                    distSqrProc.first() = Foam::sqr(GREAT);
                    distSqrProc.second() = Pstream::myProcNo();
                }
                else
                {
                    nearInfo.setPoint(treeData.centre(nearInfo.index()));

                    distSqrProc.first() = sample.distSqr(nearInfo.point());
                    distSqrProc.second() = Pstream::myProcNo();
                }
            }


            // 2. Reduce on master. Select nearest processor.

            // Find nearest - globally consistent
            Pstream::listCombineReduce(nearest, mappedPatchBase::nearestEqOp());

            // 3. Pick up my local faces that have won

            forAll(nearest, sampleI)
            {
                if (nearest[sampleI].first().hit())
                {
                    label procI = nearest[sampleI].second().second();
                    label index = nearest[sampleI].first().index();

                    if (procI == Pstream::myProcNo())
                    {
                        newPatchFaces.append(pp.addressing()[index]);
                    }
                }
            }
        }

        if (debug)
        {
            Pout<< "Found " << newPatchFaces.size()
                << " out of " << nAvailable
                << " on local processor" << endl;
        }

        patchFaces.transfer(newPatchFaces);
    }


    // Shuffle and truncate if in random mode
    const label totalSize = returnReduce(patchFaces.size(), sumOp<label>());

    if (totalSize > nAvailable)
    {
        // Check what fraction of maxPoints_ I need to generate locally.
        label myMaxPoints = scalar(patchFaces.size())/totalSize*nAvailable;

        labelList subset = identity(patchFaces.size());
        for (label iter = 0; iter < 4; ++iter)
        {
            forAll(subset, i)
            {
                label j = rndGen.position<label>(0, subset.size()-1);
                std::swap(subset[i], subset[j]);
            }
        }

        // Truncate
        subset.setSize(myMaxPoints);

        // Subset patchFaces

        if (debug)
        {
            Pout<< "In random mode : selected " << subset.size()
                << " faces out of " << patchFaces.size() << endl;
        }

        patchFaces = labelUIndList(patchFaces, subset)();
    }


    // Get points on patchFaces.
    globalIndex globalSampleNumbers(patchFaces.size());

    samplingPts.setCapacity(patchFaces.size());
    samplingCells.setCapacity(patchFaces.size());
    samplingFaces.setCapacity(patchFaces.size());
    samplingSegments.setCapacity(patchFaces.size());
    samplingCurveDist.setCapacity(patchFaces.size());

    // For calculation of min-decomp tet base points
    (void)mesh().tetBasePtIs();

    forAll(patchFaces, i)
    {
        const label facei = patchFaces[i];

        // Slightly shift point in since on warped face face-diagonal
        // decomposition might be outside cell for face-centre decomposition!
        pointIndexHit info = mappedPatchBase::facePoint
        (
            mesh(),
            facei,
            polyMesh::FACE_DIAG_TRIS
        );
        const label celli = mesh().faceOwner()[facei];

        if (info.hit())
        {
            // Move the point into the cell
            const point& cc = mesh().cellCentres()[celli];
            samplingPts.append
            (
                info.point() + 1e-1*(cc-info.point())
            );
        }
        else
        {
            samplingPts.append(info.point());
        }
        samplingCells.append(celli);
        samplingFaces.append(facei);
        samplingSegments.append(0);
        samplingCurveDist.append(globalSampleNumbers.toGlobal(i));
    }
}


void Foam::patchSeedSet::genSamples()
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

    // Move into *this
    setSamples
    (
        std::move(samplingPts),
        std::move(samplingCells),
        std::move(samplingFaces),
        std::move(samplingSegments),
        std::move(samplingCurveDist)
    );

    if (debug)
    {
        write(Info);
    }
}


void Foam::patchSeedSet::calcSamples
(
    DynamicList<point>& samplingPts,
    DynamicList<label>& samplingCells,
    DynamicList<label>& samplingFaces,
    DynamicList<label>& samplingSegments,
    DynamicList<scalar>& samplingCurveDist
)
{
    DebugInfo << "patchSeedSet : sampling on patches :" << endl;

    if (!rndGenPtr_)
    {
        rndGenPtr_.reset(new Random(0));
    }

    label nPatchPoints = 0;
    for (const label patchi : patchSet_)
    {
        const polyPatch& pp = mesh().boundaryMesh()[patchi];
        nPatchPoints += pp.size();

        DebugInfo << "    " << pp.name() << " size " << pp.size() << endl;
    }

    label nAvailable = min(maxPoints_, selectedLocations_.size());

    calcSelectedLocations
    (
        nAvailable,
        nPatchPoints,
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );

    nAvailable = maxPoints_ - nAvailable;

    calcPatchSamples
    (
        nAvailable,
        nPatchPoints,
        samplingPts,
        samplingCells,
        samplingFaces,
        samplingSegments,
        samplingCurveDist
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchSeedSet::patchSeedSet
(
    const word& name,
    const polyMesh& mesh,
    const meshSearch& searchEngine,
    const dictionary& dict
)
:
    sampledSet(name, mesh, searchEngine, dict),
    patchSet_
    (
        mesh.boundaryMesh().patchSet(dict.get<wordRes>("patches"))
    ),
    maxPoints_(dict.get<label>("maxPoints")),
    selectedLocations_
    (
        dict.getOrDefault<pointField>
        (
            "points",
            pointField(0)
        )
    )
{
    genSamples();
}


// ************************************************************************* //
