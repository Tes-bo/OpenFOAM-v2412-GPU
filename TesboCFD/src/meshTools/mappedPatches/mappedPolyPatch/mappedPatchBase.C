/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "mappedPatchBase.H"
#include "addToRunTimeSelectionTable.H"
#include "ListListOps.H"
#include "meshSearchMeshObject.H"
#include "multiWorldConnectionsObject.H"
#include "meshTools.H"
#include "OFstream.H"
#include "Random.H"
#include "treeDataFace.H"
#include "treeDataPoint.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "Time.H"
#include "mapDistribute.H"
#include "SubField.H"
#include "triangle.H"
#include "syncTools.H"
#include "treeDataCell.H"
#include "DynamicField.H"
#include "faceAreaWeightAMI.H"
#include "OTstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mappedPatchBase, 0);
}


const Foam::Enum
<
    Foam::mappedPatchBase::sampleMode
>
Foam::mappedPatchBase::sampleModeNames_
({
    { sampleMode::NEARESTCELL, "nearestCell" },
    { sampleMode::NEARESTPATCHFACE, "nearestPatchFace" },
    { sampleMode::NEARESTPATCHFACEAMI, "nearestPatchFaceAMI" },
    { sampleMode::NEARESTPATCHPOINT, "nearestPatchPoint" },
    { sampleMode::NEARESTFACE, "nearestFace" },
    { sampleMode::NEARESTONLYCELL, "nearestOnlyCell" },
});


const Foam::Enum
<
    Foam::mappedPatchBase::offsetMode
>
Foam::mappedPatchBase::offsetModeNames_
({
    { offsetMode::UNIFORM, "uniform" },
    { offsetMode::NONUNIFORM, "nonuniform" },
    { offsetMode::NORMAL, "normal" },
});


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mappedPatchBase::updateMesh(PstreamBuffers& pBufs)
{
    clearOut();
}


Foam::autoPtr<Foam::fileName> Foam::mappedPatchBase::readDatabase
(
    const dictionary& dict
)
{
    if (dict.found("sampleDatabase"))
    {
        if (dict.get<bool>("sampleDatabase"))
        {
            return autoPtr<fileName>::New
            (
                dict.getOrDefault<fileName>
                (
                    "sampleDatabasePath",
                    fileName::null
                )
            );
        }
    }
    else if (dict.found("sampleDatabasePath"))
    {
        return autoPtr<fileName>::New(dict.get<fileName>("sampleDatabasePath"));
    }

    return nullptr;
}


bool Foam::mappedPatchBase::addWorldConnection()
{
    if (sameWorld())
    {
        return true;
    }

    const Time& runTime = patch_.boundaryMesh().mesh().time();
    return const_cast<multiWorldConnections&>
    (
        multiWorldConnections::New(runTime)
    ).addConnectionByName(sampleWorld_);
}


Foam::label Foam::mappedPatchBase::getWorldCommunicator() const
{
    if (sameWorld())
    {
        return UPstream::commWorld();
    }

    const Time& runTime = patch_.boundaryMesh().mesh().time();
    return
        multiWorldConnections::New(runTime).getCommByName(sampleWorld_);
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::facePoints
(
    const polyPatch& pp
) const
{
    const polyMesh& mesh = pp.boundaryMesh().mesh();

    // Force construction of min-tet decomp
    (void)mesh.tetBasePtIs();

    // Initialise to face-centre
    auto tfacePoints = tmp<pointField>::New(patch_.size());
    auto& facePoints = tfacePoints.ref();

    forAll(pp, facei)
    {
        facePoints[facei] = facePoint
        (
            mesh,
            pp.start()+facei,
            polyMesh::FACE_DIAG_TRIS
        ).point();
    }

    return tfacePoints;
}


void Foam::mappedPatchBase::collectSamples
(
    const label mySampleWorld,      // Wanted world
    const pointField& facePoints,
    pointField& samples,            // All samples
    labelList& patchFaceWorlds,     // Per sample: wanted world
    labelList& patchFaceProcs,      // Per sample: originating processor
    labelList& patchFaces,          // Per sample: originating patchFace index
    pointField& patchFc             // Per sample: originating centre
) const
{
    DebugInFunction << nl;

    const label myComm = getCommunicator();  // Get or create
    const label myRank = Pstream::myProcNo(myComm);
    const label nProcs = Pstream::nProcs(myComm);

    const label oldWarnComm = UPstream::commWarn(myComm);

    if (debug & 2)
    {
        Perr<< "patch: " << patch_.name()
            << "[rank=" << myRank << " procs=" << nProcs
            << " comm=" << myComm << "] collect samples" << endl;
    }

    // Collect all sample points and the faces they come from.
    {
        List<pointField> globalFc(nProcs);
        globalFc[myRank] = facePoints;

        Pstream::allGatherList(globalFc, Pstream::msgType(), myComm);

        // Rework into straight list
        patchFc = ListListOps::combine<pointField>
        (
            globalFc,
            accessOp<pointField>()
        );
    }

    {
        List<pointField> globalSamples(nProcs);
        globalSamples[myRank] = samplePoints(facePoints);
        Pstream::allGatherList(globalSamples, Pstream::msgType(), myComm);

        // Rework into straight list
        samples = ListListOps::combine<pointField>
        (
            globalSamples,
            accessOp<pointField>()
        );
    }

    {
        labelListList globalFaces(nProcs);
        globalFaces[myRank] = identity(patch_.size());
        // Distribute to all processors
        Pstream::allGatherList(globalFaces, Pstream::msgType(), myComm);

        patchFaces = ListListOps::combine<labelList>
        (
            globalFaces,
            accessOp<labelList>()
        );
    }

    {
        labelList procToWorldIndex
        (
            UPstream::allGatherValues<label>(mySampleWorld, myComm)
        );
        labelList nPerProc
        (
            UPstream::allGatherValues<label>(patch_.size(), myComm)
        );

        patchFaceWorlds.resize(patchFaces.size());
        patchFaceProcs.resize(patchFaces.size());

        label sampleI = 0;
        forAll(nPerProc, proci)
        {
            for (label i = 0; i < nPerProc[proci]; i++)
            {
                patchFaceWorlds[sampleI] = procToWorldIndex[proci];
                patchFaceProcs[sampleI] = proci;
                sampleI++;
            }
        }
    }

    // Restore communicator settings
    UPstream::commWarn(oldWarnComm);
}


void Foam::mappedPatchBase::findLocalSamples
(
    const sampleMode mode,

    const label mySampleWorld,  // local world to sample == my own world
    const word& sampleRegion,   // local region to sample
    const word& samplePatch,    // local patch to sample

    const pointField& samples,
    List<nearInfoWorld>& nearest
) const
{
    DebugInFunction << nl;

    const label myComm = getCommunicator();  // Get or create

    // Find the local cell containing the samples
    const label myRank = Pstream::myProcNo(myComm);

    // Lookup the correct region
    const polyMesh& mesh = lookupMesh(sampleRegion);

    // All the info for nearest. Construct to miss
    nearest.setSize(samples.size());
    nearInfoWorld miss;
    {
        miss.first().second() = Tuple2<scalar, label>(Foam::sqr(GREAT), -1);
        miss.second() = -1; // set world to be ignored
    }
    nearest = miss;

    switch (mode)
    {
        case NEARESTCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                label celli = tree.findInside(sample);

                if (celli == -1)
                {
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
                else
                {
                    const point& cc = mesh.cellCentres()[celli];

                    near.first().first() = pointIndexHit
                    (
                        true,
                        cc,
                        celli
                    );
                    near.first().second().first() = sample.distSqr(cc);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            break;
        }

        case NEARESTONLYCELL:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const indexedOctree<Foam::treeDataCell>& tree = mesh.cellTree();

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                near.first().first() = tree.findNearest(sample, sqr(GREAT));

                near.first().second().first() =
                    near.first().first().hitPoint().distSqr(sample);

                near.first().second().second() = myRank;
                near.second() = mySampleWorld;
            }
            break;
        }

        case NEARESTPATCHFACE:
        {
            Random rndGen(123456);

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearInfoWorld& near = nearest[sampleI];
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            else
            {
                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1e-4,
                        ROOTVSMALL
                    )
                );

                indexedOctree<treeDataFace> boundaryTree
                (
                    treeDataFace(mesh, pp.range()),  // Patch faces

                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );
                const auto& treeData = boundaryTree.shapes();

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    nearInfoWorld& near = nearest[sampleI];
                    pointIndexHit& nearInfo = near.first().first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        patchBb.magSqr()
                    );

                    if (!nearInfo.hit())
                    {
                        near.first().second().first() = Foam::sqr(GREAT);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                    else
                    {
                        const point& fc = treeData.centre(nearInfo.index());

                        nearInfo.setPoint(fc);
                        near.first().second().first() = sample.distSqr(fc);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                }
            }
            break;
        }

        case NEARESTPATCHPOINT:
        {
            Random rndGen(123456);

            const polyPatch& pp = lookupPatch(sampleRegion, samplePatch);

            if (pp.empty())
            {
                forAll(samples, sampleI)
                {
                    nearInfoWorld& near = nearest[sampleI];
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            else
            {
                // patch (local) points
                treeBoundBox patchBb
                (
                    treeBoundBox(pp.points(), pp.meshPoints()).extend
                    (
                        rndGen,
                        1e-4,
                        ROOTVSMALL
                    )
                );

                indexedOctree<treeDataPoint> boundaryTree
                (
                    // Patch points
                    treeDataPoint(mesh.points(), pp.meshPoints()),

                    patchBb,        // overall search domain
                    8,              // maxLevel
                    10,             // leafsize
                    3.0             // duplicity
                );

                forAll(samples, sampleI)
                {
                    const point& sample = samples[sampleI];

                    nearInfoWorld& near = nearest[sampleI];
                    pointIndexHit& nearInfo = near.first().first();
                    nearInfo = boundaryTree.findNearest
                    (
                        sample,
                        patchBb.magSqr()
                    );

                    if (!nearInfo.hit())
                    {
                        near.first().second().first() = Foam::sqr(GREAT);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                    else
                    {
                        const point& pt = nearInfo.point();

                        near.first().second().first() = sample.distSqr(pt);
                        near.first().second().second() = myRank;
                        near.second() = mySampleWorld;
                    }
                }
            }
            break;
        }

        case NEARESTFACE:
        {
            if (samplePatch.size() && samplePatch != "none")
            {
                FatalErrorInFunction
                    << "No need to supply a patch name when in "
                    << sampleModeNames_[mode] << " mode." << exit(FatalError);
            }

            //- Note: face-diagonal decomposition
            const meshSearchMeshObject& meshSearchEngine =
                meshSearchMeshObject::New(mesh);

            forAll(samples, sampleI)
            {
                const point& sample = samples[sampleI];
                nearInfoWorld& near = nearest[sampleI];

                label facei = meshSearchEngine.findNearestFace(sample);

                if (facei == -1)
                {
                    near.first().second().first() = Foam::sqr(GREAT);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
                else
                {
                    const point& fc = mesh.faceCentres()[facei];

                    near.first().first() = pointIndexHit(true, fc, facei);
                    near.first().second().first() = sample.distSqr(fc);
                    near.first().second().second() = myRank;
                    near.second() = mySampleWorld;
                }
            }
            break;
        }

        case NEARESTPATCHFACEAMI:
        {
            // nothing to do here
            return;
        }

        default:
        {
            FatalErrorInFunction
                << "problem." << abort(FatalError);
        }
    }
}


void Foam::mappedPatchBase::findSamples
(
    const sampleMode mode,
    const label myWorld,
    const pointField& samples,
    const labelList& wantedWorlds,
    const labelList& origProcs,

    labelList& sampleProcs,
    labelList& sampleIndices,
    pointField& sampleLocations
) const
{
    DebugInFunction << nl;

    // Find the processor/cell containing the samples. Does not account
    // for samples being found in two processors.

    const label myComm = getCommunicator();  // Get or create
    const label myRank = Pstream::myProcNo(myComm);
    const label nProcs = Pstream::nProcs(myComm);

    const label oldWarnComm = UPstream::commWarn(myComm);

    wordList samplePatches(nProcs);
    {
        samplePatches[myRank] = samplePatch_;
        Pstream::allGatherList(samplePatches, Pstream::msgType(), myComm);
    }
    wordList sampleRegions(nProcs);
    {
        sampleRegions[myRank] = sampleRegion_;
        Pstream::allGatherList(sampleRegions, Pstream::msgType(), myComm);
    }


    // Find all the info for nearest
    List<nearInfoWorld> nearest(samples.size());
    forAll(nearest, samplei)
    {
        nearest[samplei].first() = nearInfo
        (
            pointIndexHit(),
            Tuple2<scalar, label>(Foam::sqr(GREAT), -1)
        );
        nearest[samplei].second() = wantedWorlds[samplei];
    }


    // Extract samples to search for locally
    {
        DynamicList<label> localMap(samples.size());
        forAll(wantedWorlds, samplei)
        {
            if (wantedWorlds[samplei] == myWorld)
            {
                localMap.append(samplei);
            }
        }

        if (localMap.size())
        {
            pointField localSamples(samples, localMap);
            labelList localOrigProcs(origProcs, localMap);

            //Assume single patch to sample for now
            const word localOrigPatch(samplePatches[localOrigProcs[0]]);
            const word localOrigRegion(sampleRegions[localOrigProcs[0]]);
            List<nearInfoWorld> localNearest(localSamples.size());

            if (debug)
            {
                Pout<< "Searching locally for " << localSamples.size()
                    << " samples on region:" << localOrigRegion
                    << " on patch:" << localOrigPatch << endl;
            }
            findLocalSamples
            (
                mode,
                myWorld,
                localOrigRegion,
                localOrigPatch,
                localSamples,
                localNearest
            );
            UIndirectList<nearInfoWorld>(nearest, localMap) = localNearest;
        }
    }


    // Find nearest - globally consistent
    Pstream::listCombineReduce
    (
        nearest,
        nearestWorldEqOp(),
        Pstream::msgType(),
        myComm
    );

    //if (debug)
    //{
    //    Pout<< "** After combining:" << endl;
    //    forAll(nearest, samplei)
    //    {
    //        Pout<< "  sample:" << samples[samplei]
    //            << " origating from proc:" << origProcs[samplei]
    //            << " to be found on world:" << wantedWorlds[samplei] << nl
    //            << "    found on world:" << nearest[samplei].second() << nl
    //            << "    found on proc:"
    //            << nearest[samplei].first().second().second() << nl
    //            << "    found on patchfacei:"
    //            << nearest[samplei].first().first().index() << nl
    //            << "    found at location:"
    //            << nearest[samplei].first().first().point() << nl;
    //    }
    //    Pout<< endl;
    //}

    // Convert back into proc+local index
    sampleProcs.setSize(samples.size());
    sampleIndices.setSize(samples.size());
    sampleLocations.setSize(samples.size());

    forAll(nearest, sampleI)
    {
        const nearInfo& ni = nearest[sampleI].first();

        if (!ni.first().hit())
        {
            sampleProcs[sampleI] = -1;
            sampleIndices[sampleI] = -1;
            sampleLocations[sampleI] = vector::max;
        }
        else
        {
            sampleProcs[sampleI] = ni.second().second();
            sampleIndices[sampleI] = ni.first().index();
            sampleLocations[sampleI] = ni.first().point();
        }
    }

    // Return communicator settings
    UPstream::commWarn(oldWarnComm);
}


void Foam::mappedPatchBase::calcMapping() const
{
    DebugInFunction;

    static bool hasWarned = false;
    if (mapPtr_)
    {
        FatalErrorInFunction
            << "Mapping already calculated" << exit(FatalError);
    }

    DebugInFunction << nl;


    // Early construction of tetBasePtIs since does parallel
    // communication which might conflict with inter-world communicator
    // (for multi-world operation)
    (void)patch_.boundaryMesh().mesh().tetBasePtIs();

    const label myComm = getCommunicator();  // Get or create

    //// Make sure if running in database that there is a syncObjects FO
    //if (sampleDatabase() && !sameWorld())
    //{
    //    const word syncName("syncObjects");
    //    const polyMesh& mesh = patch_.boundaryMesh().mesh();
    //    const Time& runTime = mesh.time();
    //    const functionObjectList& fos = runTime.functionObjects();
    //
    //    forAll(fos, i)
    //    {
    //        Pout<< "** FO:" << fos[i].name() << " tpye:" << fos[i].type()
    //            << endl;
    //    }
    //
    //    if (fos.findObjectID(syncName) == -1)
    //    {
    //        //FatalErrorInFunction
    //        WarningInFunction
    //            << "Did not detect functionObject " << syncName
    //            << ". This is used to synchronise data inbetween worlds"
    //            << endl;
    //    }
    //}


    // Get points on face (since cannot use face-centres - might be off
    // face-diagonal decomposed tets.
    tmp<pointField> patchPoints(facePoints(patch_));

    // Get offsetted points
    const pointField offsettedPoints(samplePoints(patchPoints()));

    // Do a sanity check - am I sampling my own patch?
    // This only makes sense for a non-zero offset.
    bool sampleMyself =
    (
        mode_ == NEARESTPATCHFACE
     && sampleWorld() == UPstream::myWorld()
     && sampleRegion() == patch_.boundaryMesh().mesh().name()
     && samplePatch() == patch_.name()
    );

    if (sampleMyself)
    {
        // Check offset
        vectorField d(offsettedPoints-patchPoints());
        bool coincident = (gAverage(mag(d)) <= ROOTVSMALL);

        if (sampleMyself && coincident)
        {
            WarningInFunction
                << "Invalid offset " << d << endl
                << "Offset is the vector added to the patch face centres to"
                << " find the patch face supplying the data." << endl
                << "Setting it to " << d
                << " on the same patch, on the same region"
                << " will find the faces themselves which does not make sense"
                << " for anything but testing." << endl
                << "patch:" << patch_.name() << endl
                << "sampleRegion:" << sampleRegion() << endl
                << "mode:" << sampleModeNames_[mode_] << endl
                << "samplePatch:" << samplePatch() << endl
                << "offsetMode:" << offsetModeNames_[offsetMode_] << endl;
        }
    }

    // Get local world and world-to-sample in index form
    const label myWorld = Pstream::myWorldID();
    const label mySampleWorld = Pstream::allWorlds().find(sampleWorld_);


    // Get global list of all samples and the processor and face they come from.
    pointField samples;         // coordinates
    labelList patchFaceWorlds;  // world to sample
    labelList patchFaceProcs;   // originating processor
    labelList patchFaces;       // originating face on processor patch
    pointField patchFc;         // originating face centre
    collectSamples
    (
        mySampleWorld,          // world I want to sample
        patchPoints,

        samples,
        patchFaceWorlds,
        patchFaceProcs,
        patchFaces,
        patchFc
    );

    //if (debug)
    //{
    //    forAll(samples, samplei)
    //    {
    //        Pout<< "    sample:" << samples[samplei]
    //            << " origating from proc:" << patchFaceProcs[samplei]
    //            << "  face:" << patchFaces[samplei]
    //            << " to be found on world:" << patchFaceWorlds[samplei] << nl;
    //    }
    //}

    // Find processor and cell/face samples are in and actual location.
    labelList sampleProcs;
    labelList sampleIndices;
    pointField sampleLocations;
    findSamples
    (
        mode_,
        myWorld,        // my world (in index form)
        samples,
        patchFaceWorlds,
        patchFaceProcs,

        sampleProcs,
        sampleIndices,
        sampleLocations
    );

    // Check for samples that were not found. This will only happen for
    // NEARESTCELL since finds cell containing a location
    if (mode_ == NEARESTCELL)
    {
        label nNotFound = 0;
        forAll(sampleProcs, sampleI)
        {
            if (sampleProcs[sampleI] == -1)
            {
                nNotFound++;
            }
        }
        reduce(nNotFound, sumOp<label>(), Pstream::msgType(), myComm);

        if (nNotFound > 0)
        {
            if (!hasWarned)
            {
                WarningInFunction
                    << "Did not find " << nNotFound
                    << " out of " << sampleProcs.size() << " total samples."
                    << " Sampling these on owner cell centre instead." << endl
                    << "On patch " << patch_.name()
                    << " on region " << sampleRegion()
                    << " in mode " << sampleModeNames_[mode_] << endl
                    << "with offset mode " << offsetModeNames_[offsetMode_]
                    << ". Suppressing further warnings from " << type() << endl;

                hasWarned = true;
            }

            // Collect the samples that cannot be found
            DynamicList<label> subMap;
            forAll(sampleProcs, sampleI)
            {
                if (sampleProcs[sampleI] == -1)
                {
                    subMap.append(sampleI);
                }
            }

            // And re-search for pure nearest (should not fail)
            labelList subSampleProcs;
            labelList subSampleIndices;
            pointField subSampleLocations;
            findSamples
            (
                NEARESTONLYCELL,
                myWorld,        // my world (in index form)

                pointField(samples, subMap),
                UIndirectList<label>(patchFaceWorlds, subMap)(),
                UIndirectList<label>(patchFaceProcs, subMap)(),

                subSampleProcs,
                subSampleIndices,
                subSampleLocations
            );

            // Insert
            labelUIndList(sampleProcs, subMap) = subSampleProcs;
            labelUIndList(sampleIndices, subMap) = subSampleIndices;
            UIndirectList<point>(sampleLocations, subMap) = subSampleLocations;
        }
    }

    // Now we have all the data we need:
    // - where sample originates from (so destination when mapping):
    //   patchFaces, patchFaceProcs.
    // - cell/face sample is in (so source when mapping)
    //   sampleIndices, sampleProcs.

    if (Pstream::master(myComm))
    {
        forAll(samples, i)
        {
            if (sampleProcs[i] == -1)
            {
                FatalErrorInFunction << "did not find sample "
                    << patchFc[i] << " on patch " << patch_.name()
                    << " on region "
                    << patch_.boundaryMesh().mesh().name()
                    << " on processor " << patchFaceProcs[i]
                    << exit(FatalError);
            }
        }
    }


    if (debug && Pstream::master(myComm))
    {
        //forAll(samples, i)
        //{
        //    Pout<< i << " need data in region "
        //        << patch_.boundaryMesh().mesh().name()
        //        << " for proc:" << patchFaceProcs[i]
        //        << " face:" << patchFaces[i]
        //        << " at:" << patchFc[i] << endl
        //        << "Found data in region " << sampleRegion()
        //        << " at proc:" << sampleProcs[i]
        //        << " face:" << sampleIndices[i]
        //        << " at:" << sampleLocations[i]
        //        << nl << endl;
        //}

        OFstream str
        (
            patch_.boundaryMesh().mesh().time().path()
          / patch_.name()
          + "_mapped.obj"
        );
        Pout<< "Dumping mapping as lines from patch faceCentres to"
            << " sampled cell/faceCentres/points to file " << str.name()
            << endl;

        label vertI = 0;

        forAll(patchFc, i)
        {
            meshTools::writeOBJ(str, patchFc[i]);
            vertI++;
            meshTools::writeOBJ(str, sampleLocations[i]);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }

    // Determine schedule.
    mapPtr_.reset(new mapDistribute(sampleProcs, patchFaceProcs, myComm));

    // Rework the schedule from indices into samples to cell data to send,
    // face data to receive.

    labelListList& subMap = mapPtr_().subMap();
    labelListList& constructMap = mapPtr_().constructMap();

    forAll(subMap, proci)
    {
        subMap[proci] = labelUIndList(sampleIndices, subMap[proci]);
        constructMap[proci] = labelUIndList(patchFaces, constructMap[proci]);

        if (debug)
        {
            Pout<< "To proc:" << proci << " sending values of cells/faces:"
                << subMap[proci] << endl;
            Pout<< "From proc:" << proci
                << " receiving values of patch faces:"
                << constructMap[proci] << endl;
        }
    }


    // Redo constructSize
    mapPtr_().constructSize() = patch_.size();

    if (debug)
    {
        // Check that all elements get a value.
        bitSet used(patch_.size());
        forAll(constructMap, proci)
        {
            const labelList& map = constructMap[proci];

            forAll(map, i)
            {
                label facei = map[i];

                if (used.test(facei))
                {
                    FatalErrorInFunction
                        << "On patch " << patch_.name()
                        << " patchface " << facei
                        << " is assigned to more than once."
                        << abort(FatalError);
                }
                else
                {
                    used.set(facei);
                }
            }
        }
        forAll(used, facei)
        {
            if (!used.test(facei))
            {
                FatalErrorInFunction
                    << "On patch " << patch_.name()
                    << " patchface " << facei
                    << " is never assigned to."
                    << abort(FatalError);
            }
        }
    }

    updateMeshTime().setUpToDate();
    if (sameWorld())
    {
        updateSampleMeshTime().setUpToDate();
    }
}


const Foam::autoPtr<Foam::searchableSurface>& Foam::mappedPatchBase::surfPtr()
const
{
    const word surfType(surfDict_.getOrDefault<word>("type", "none"));

    if (!surfPtr_ && surfType != "none")
    {
        word surfName(surfDict_.getOrDefault("name", patch_.name()));

        const polyMesh& mesh = patch_.boundaryMesh().mesh();

        surfPtr_ =
            searchableSurface::New
            (
                surfType,
                IOobject
                (
                    surfName,
                    mesh.time().constant(),
                    "triSurface",
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                surfDict_
            );
    }

    return surfPtr_;
}


void Foam::mappedPatchBase::calcAMI() const
{
    if (AMIPtr_->upToDate())
    {
        DebugInFunction
            << "AMI already up-to-date"
            << endl;

        return;
    }

    DebugInFunction << nl;

    const label myComm = getCommunicator();  // Get or create

    // Pre-calculate surface (if any)
    const auto& surf = surfPtr();

    // Check if running locally
    if (sampleWorld_.empty() || sameWorld())
    {
        const polyPatch& nbr = samplePolyPatch();

        // Transform neighbour patch to local system
        const pointField nbrPoints(samplePoints(nbr.localPoints()));

        const primitivePatch nbrPatch0
        (
            SubList<face>
            (
                nbr.localFaces(),
                nbr.size()
            ),
            nbrPoints
        );


        if (debug)
        {
            OFstream os(patch_.name() + "_neighbourPatch-org.obj");
            meshTools::writeOBJ(os, samplePolyPatch().localFaces(), nbrPoints);

            OFstream osN(patch_.name() + "_neighbourPatch-trans.obj");
            meshTools::writeOBJ(osN, nbrPatch0, nbrPoints);

            OFstream osO(patch_.name() + "_ownerPatch.obj");
            meshTools::writeOBJ(osO, patch_.localFaces(), patch_.localPoints());
        }

        // Construct/apply AMI interpolation to determine addressing and
        // weights.

        // Change to use inter-world communicator
        const label oldWarnComm = UPstream::commWarn(myComm);
        const label oldWorldComm = UPstream::commWorld(myComm);

        AMIPtr_->calculate(patch_, nbrPatch0, surf);

        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);
    }
    else
    {
        faceList dummyFaces;
        pointField dummyPoints;
        const primitivePatch dummyPatch
        (
            SubList<face>(dummyFaces),
            dummyPoints
        );

        // Change to use inter-world communicator
        AMIPtr_->comm(myComm);

        const label oldWarnComm = UPstream::commWarn(myComm);
        const label oldWorldComm = UPstream::commWorld(myComm);

        if (masterWorld())
        {
            // Construct/apply AMI interpolation to determine addressing
            // and weights. Have patch_ for src faces, 0 faces for the
            // target side
            AMIPtr_->calculate(patch_, dummyPatch, surf);
        }
        else
        {
            // Construct/apply AMI interpolation to determine addressing
            // and weights. Have 0 faces for src side, patch_ for the tgt
            // side
            AMIPtr_->calculate(dummyPatch, patch_, surf);
        }
        // Now the AMI addressing/weights will be from src side (on masterWorld
        // processors) to tgt side (on other processors)

        // Restore communicator settings
        UPstream::commWarn(oldWarnComm);
        UPstream::commWorld(oldWorldComm);
    }
}


const Foam::objectRegistry& Foam::mappedPatchBase::subRegistry
(
    const objectRegistry& obr,
    const wordList& names,
    const label index
)
{
    const objectRegistry& sub = obr.subRegistry(names[index], true);
    if (index == names.size()-1)
    {
        return sub;
    }
    else
    {
        return subRegistry(sub, names, index+1);
    }
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::mappedPatchBase(const polyPatch& pp)
:
    patch_(pp),
    sampleWorld_(),
    sampleRegion_(patch_.boundaryMesh().mesh().name()),
    mode_(NEARESTPATCHFACE),
    samplePatch_(),
    coupleGroup_(),
    sampleDatabasePtr_(),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(pp.size(), offset_),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_(true),
    mapPtr_(nullptr),
    AMIReverse_(false),
    AMIPtr_(new faceAreaWeightAMI(true, AMIReverse_)),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    // NOTE: same region, no sample-world. Thus no world-world communication
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vectorField& offsets
)
:
    mappedPatchBase
    (
        pp,
        sampleRegion,
        mode,
        samplePatch,
        scalar(0)
    )
{
    mappedPatchBase::setOffset(offsets);
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const vector& uniformOffset
)
:
    mappedPatchBase
    (
        pp,
        sampleRegion,
        mode,
        samplePatch,
        scalar(0)
    )
{
    mappedPatchBase::setOffset(uniformOffset);
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const word& sampleRegion,
    const sampleMode mode,
    const word& samplePatch,
    const scalar normalDistance
)
:
    patch_(pp),
    sampleWorld_(),
    sampleRegion_(sampleRegion),
    mode_(mode),
    samplePatch_(samplePatch),
    coupleGroup_(),
    sampleDatabasePtr_(),
    offsetMode_(NORMAL),
    offset_(Zero),
    offsets_(0),
    distance_(normalDistance),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_(false),
    AMIPtr_(new faceAreaWeightAMI(true, AMIReverse_)),
    surfPtr_(nullptr),
    surfDict_(fileName("surface"))
{
    addWorldConnection();
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.getOrDefault<word>("sampleWorld", word::null)),
    sampleRegion_(dict.getOrDefault<word>("sampleRegion", word::null)),
    mode_(sampleModeNames_.get("sampleMode", dict)),
    samplePatch_(dict.getOrDefault<word>("samplePatch", word::null)),
    coupleGroup_(dict),
    sampleDatabasePtr_(readDatabase(dict)),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_
    (
        dict.getOrDefault
        (
            "reverseTarget",        // AMIInterpolation uses this keyword
            dict.getOrDefault
            (
                "flipNormals",
                false
            )
        )
    ),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            AMIReverse_
        )
    ),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    addWorldConnection();

    if (!coupleGroup_.good())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }

    if (offsetModeNames_.readIfPresent("offsetMode", dict, offsetMode_))
    {
        switch (offsetMode_)
        {
            case UNIFORM:
            {
                dict.readEntry("offset", offset_);
            }
            break;

            case NONUNIFORM:
            {
                offsets_ = pointField("offsets", dict, patch_.size());
            }
            break;

            case NORMAL:
            {
                dict.readEntry("distance", distance_);
            }
            break;
        }
    }
    else if (dict.readIfPresent("offset", offset_))
    {
        offsetMode_ = UNIFORM;
    }
    else if (dict.found("offsets"))
    {
        offsetMode_ = NONUNIFORM;
        offsets_ = pointField("offsets", dict, patch_.size());
    }
    else if (mode_ != NEARESTPATCHFACE && mode_ != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Please supply the offsetMode as one of "
            << offsetModeNames_
            << exit(FatalIOError);
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const sampleMode mode,
    const dictionary& dict
)
:
    patch_(pp),
    sampleWorld_(dict.getOrDefault<word>("sampleWorld", word::null)),
    sampleRegion_(dict.getOrDefault<word>("sampleRegion", word::null)),
    mode_(mode),
    samplePatch_(dict.getOrDefault<word>("samplePatch", word::null)),
    coupleGroup_(dict), //dict.getOrDefault<word>("coupleGroup", word::null)),
    sampleDatabasePtr_(readDatabase(dict)),
    offsetMode_(UNIFORM),
    offset_(Zero),
    offsets_(0),
    distance_(0),
    communicator_(-1),  // Demand-driven (cached value)
    sameRegion_
    (
        sampleWorld_.empty()
     && sampleRegion_ == patch_.boundaryMesh().mesh().name()
    ),
    mapPtr_(nullptr),
    AMIReverse_
    (
        dict.getOrDefault
        (
            "reverseTarget",        // AMIInterpolation uses this keyword
            dict.getOrDefault
            (
                "flipNormals",
                false
            )
        )
    ),
    AMIPtr_
    (
        AMIInterpolation::New
        (
            dict.getOrDefault("AMIMethod", faceAreaWeightAMI::typeName),
            dict,
            AMIReverse_
        )
    ),
    surfPtr_(nullptr),
    surfDict_(dict.subOrEmptyDict("surface"))
{
    addWorldConnection();

    if (mode != NEARESTPATCHFACE && mode != NEARESTPATCHFACEAMI)
    {
        FatalIOErrorInFunction(dict)
            << "Construct from sampleMode and dictionary only applicable for "
            << " collocated patches in modes "
            << sampleModeNames_[NEARESTPATCHFACE] << ','
            << sampleModeNames_[NEARESTPATCHFACEAMI]
            << exit(FatalIOError);
    }

    if (!coupleGroup_.good())
    {
        if (sampleWorld_.empty() && sampleRegion_.empty())
        {
            // If no coupleGroup and no sampleRegion assume local region
            sampleRegion_ = patch_.boundaryMesh().mesh().name();
            sameRegion_ = true;
        }
    }
}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    sampleDatabasePtr_
    (
        mpb.sampleDatabasePtr_
      ? new fileName(mpb.sampleDatabasePtr_())
      : nullptr
    ),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_(mpb.offsets_),
    distance_(mpb.distance_),
    communicator_(mpb.communicator_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    AMIPtr_(mpb.AMIPtr_->clone()),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


Foam::mappedPatchBase::mappedPatchBase
(
    const polyPatch& pp,
    const mappedPatchBase& mpb,
    const labelUList& mapAddressing
)
:
    patch_(pp),
    sampleWorld_(mpb.sampleWorld_),
    sampleRegion_(mpb.sampleRegion_),
    mode_(mpb.mode_),
    samplePatch_(mpb.samplePatch_),
    coupleGroup_(mpb.coupleGroup_),
    sampleDatabasePtr_
    (
        mpb.sampleDatabasePtr_
      ? new fileName(mpb.sampleDatabasePtr_())
      : nullptr
    ),
    offsetMode_(mpb.offsetMode_),
    offset_(mpb.offset_),
    offsets_
    (
        offsetMode_ == NONUNIFORM
      ? vectorField(mpb.offsets_, mapAddressing)
      : vectorField()
    ),
    distance_(mpb.distance_),
    communicator_(mpb.communicator_),
    sameRegion_(mpb.sameRegion_),
    mapPtr_(nullptr),
    AMIReverse_(mpb.AMIReverse_),
    AMIPtr_(mpb.AMIPtr_->clone()),
    surfPtr_(nullptr),
    surfDict_(mpb.surfDict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mappedPatchBase::~mappedPatchBase()
{
    clearOut();
}


void Foam::mappedPatchBase::clearOut()
{
    mapPtr_.reset(nullptr);
    AMIPtr_->upToDate(false);

    //Note: no need to clear out surface since not mesh related
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedPatchBase::setOffset(const scalar normalDist)
{
    clearOut();
    offsetMode_ = offsetMode::NORMAL;
    offset_ = Zero;
    offsets_.clear();
    distance_ = normalDist;
}


void Foam::mappedPatchBase::setOffset(const vector& uniformOffset)
{
    clearOut();
    offsetMode_ = offsetMode::UNIFORM;
    offset_ = uniformOffset;
    offsets_.clear();
    distance_ = Zero;
}


void Foam::mappedPatchBase::setOffset(const vectorField& offsets)
{
    clearOut();
    offsetMode_ = offsetMode::NONUNIFORM;
    offset_ = Zero;
    offsets_ = offsets;
    distance_ = Zero;
}


const Foam::polyMesh& Foam::mappedPatchBase::lookupMesh
(
    const word& sampleRegion
) const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();
    return
    (
        sampleRegion.empty() || sampleRegion == thisMesh.name()
      ? thisMesh
      : thisMesh.time().lookupObject<polyMesh>(sampleRegion)
    );
}


const Foam::polyPatch& Foam::mappedPatchBase::lookupPatch
(
    const word& sampleRegion,
    const word& samplePatch
) const
{
    const polyMesh& nbrMesh = lookupMesh(sampleRegion);

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch
            << " in region " << sampleRegion_ << endl
            << exit(FatalError);
    }
    return nbrMesh.boundaryMesh()[patchi];
}


const Foam::polyMesh& Foam::mappedPatchBase::sampleMesh() const
{
    if (UPstream::myWorld() != sampleWorld_)
    {
        FatalErrorInFunction
            << "sampleWorld : " << sampleWorld_
            << " is not the current world : " << UPstream::myWorld()
            << exit(FatalError);
    }
    return lookupMesh(sampleRegion());
}


const Foam::polyPatch& Foam::mappedPatchBase::samplePolyPatch() const
{
    const polyMesh& nbrMesh = sampleMesh();

    const label patchi = nbrMesh.boundaryMesh().findPatchID(samplePatch());

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << samplePatch()
            << " in region " << sampleRegion_ << endl
            << "Valid patches are " << nbrMesh.boundaryMesh().names()
            << exit(FatalError);
    }

    return nbrMesh.boundaryMesh()[patchi];
}


bool Foam::mappedPatchBase::upToDate() const
{
    const polyMesh& thisMesh = patch_.boundaryMesh().mesh();

    bool thisUpToDate = thisMesh.upToDatePoints(updateMeshTime());
    bool sampleUpToDate =
    (
        sameWorld()
      ? sampleMesh().upToDatePoints(updateSampleMeshTime())
      : true
    );


    // Lambda to check for points on the patch being the same
    auto checkPointMovement = []
    (
        const polyMesh& mesh,
        const polyPatch& patch,
        regIOobject& state
    ) -> bool
    {
        bool upToDate = true;
        const auto& oldPoints = mesh.oldPoints();
        const auto& points = mesh.points();

        for (const label pointi : patch.meshPoints())
        {
            const point& oldPt = oldPoints[pointi];
            const point& currentPt = points[pointi];

            if (mag(oldPt - currentPt) > SMALL)
            {
                upToDate = false;
                break;
            }
        }

        Pstream::reduceAnd(upToDate);

        if (upToDate)
        {
            state.setUpToDate();
        }

        return upToDate;
    };


    if (!thisUpToDate && thisMesh.moving())
    {
        // Moving (but not topoChanging mesh) : do more accurate check:
        // compare actual patch point position

        thisUpToDate = checkPointMovement(thisMesh, patch_, updateMeshTime());
    }

    if
    (
       !sampleUpToDate
     && sampleMesh().moving()
     && (
            mode_ == NEARESTPATCHFACE
         || mode_ == NEARESTPATCHFACEAMI
         || mode_ == NEARESTPATCHPOINT
        )
    )
    {
        sampleUpToDate = checkPointMovement
        (
            sampleMesh(),
            samplePolyPatch(),
            updateSampleMeshTime()
        );
    }

    return (thisUpToDate && sampleUpToDate);
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints
(
    const pointField& fc
) const
{
    auto tfld = tmp<pointField>::New(fc);
    auto& fld = tfld.ref();

    switch (offsetMode_)
    {
        case UNIFORM:
        {
            fld += offset_;
            break;
        }
        case NONUNIFORM:
        {
            fld += offsets_;
            break;
        }
        case NORMAL:
        {
            // Get outwards pointing normal
            vectorField n(patch_.faceAreas());
            n /= mag(n);

            fld += distance_*n;
            break;
        }
    }

    return tfld;
}


Foam::tmp<Foam::pointField> Foam::mappedPatchBase::samplePoints() const
{
    return samplePoints(facePoints(patch_));
}


Foam::pointIndexHit Foam::mappedPatchBase::facePoint
(
    const polyMesh& mesh,
    const label facei,
    const polyMesh::cellDecomposition decompMode
)
{
    const point& fc = mesh.faceCentres()[facei];

    switch (decompMode)
    {
        case polyMesh::FACE_PLANES:
        case polyMesh::FACE_CENTRE_TRIS:
        {
            // For both decompositions the face centre is guaranteed to be
            // on the face
            return pointIndexHit(true, fc, facei);
        }
        break;

        case polyMesh::FACE_DIAG_TRIS:
        case polyMesh::CELL_TETS:
        {
            // Find the intersection of a ray from face centre to cell centre
            // Find intersection of (face-centre-decomposition) centre to
            // cell-centre with face-diagonal-decomposition triangles.

            const pointField& p = mesh.points();
            const face& f = mesh.faces()[facei];

            if (f.size() <= 3)
            {
                // Return centre of triangle.
                return pointIndexHit(true, fc, 0);
            }

            const label celli = mesh.faceOwner()[facei];
            const point& cc = mesh.cellCentres()[celli];
            vector d = fc-cc;

            const label fp0 = mesh.tetBasePtIs()[facei];
            const point& basePoint = p[f[fp0]];

            label fp = f.fcIndex(fp0);
            for (label i = 2; i < f.size(); i++)
            {
                const point& thisPoint = p[f[fp]];
                label nextFp = f.fcIndex(fp);
                const point& nextPoint = p[f[nextFp]];

                const triPointRef tri(basePoint, thisPoint, nextPoint);
                pointHit hitInfo = tri.intersection
                (
                    cc,
                    d,
                    intersection::HALF_RAY
                );

                if (hitInfo.hit() && hitInfo.distance() > 0)
                {
                    return pointIndexHit(true, hitInfo.point(), i-2);
                }

                fp = nextFp;
            }

            // Fall-back
            return pointIndexHit(false, fc, -1);
        }
        break;

        default:
        {
            FatalErrorInFunction
                << "problem" << abort(FatalError);
            return pointIndexHit();
        }
    }
}


const Foam::objectRegistry& Foam::mappedPatchBase::subRegistry
(
    const objectRegistry& obr,
    const fileName& path
)
{
    // Lookup (and create if non-existing) a registry using
    // '/' separated path. Like 'mkdir -p'

    fileName cleanedPath(path);
    cleanedPath.clean();  // Remove unneeded ".."
    const wordList names(cleanedPath.components());

    if (names.empty())
    {
        return obr;
    }
    else
    {
        return subRegistry(obr, names, 0);
    }
}


Foam::fileName Foam::mappedPatchBase::sendPath
(
    const fileName& root,
    const label proci
)
{
    const word processorName("processor"+Foam::name(proci));
    return root/"send"/processorName;
}


Foam::fileName Foam::mappedPatchBase::sendPath(const label proci) const
{
    return sendPath(sampleDatabasePath(), proci);
}


Foam::fileName Foam::mappedPatchBase::receivePath
(
    const fileName& root,
    const label proci
)
{
    const word processorName("processor"+Foam::name(proci));
    return root/"receive"/processorName;
}


Foam::fileName Foam::mappedPatchBase::receivePath(const label proci) const
{
    return receivePath(sampleDatabasePath(), proci);
}


void Foam::mappedPatchBase::writeDict
(
    const objectRegistry& obr,
    dictionary& dict
)
{
    forAllIters(obr, iter)
    {
        regIOobject* objPtr = iter.val();
        const regIOobject& obj = *objPtr;

        if (isA<objectRegistry>(obj))
        {
            dictionary& d = dict.subDictOrAdd(obj.name());
            writeDict(dynamic_cast<const objectRegistry&>(obj), d);
        }
        else if
        (
            writeIOField<scalar>(obj, dict)
         || writeIOField<vector>(obj, dict)
         || writeIOField<sphericalTensor>(obj, dict)
         || writeIOField<symmTensor>(obj, dict)
         || writeIOField<tensor>(obj, dict)
        )
        {
            // IOField specialisation
        }
        else
        {
            // Not tested. No way of retrieving data anyway ...
            OTstream os;
            obj.writeData(os);

            primitiveEntry* pePtr = new primitiveEntry(obj.name(), os.tokens());
            dict.add(pePtr);
        }
    }
}


void Foam::mappedPatchBase::readDict(const dictionary& d, objectRegistry& obr)
{
    // Reverse of writeDict - reads fields from dictionary into objectRegistry
    for (const entry& e : d)
    {
        if (e.isDict())
        {
            // Add sub registry
            objectRegistry& sub = const_cast<objectRegistry&>
            (
                obr.subRegistry(e.keyword(), true)
            );

            readDict(e.dict(), sub);
        }
        else
        {
            ITstream& is = e.stream();
            token tok(is);

            if
            (
                constructIOField<scalar>(e.keyword(), tok, is, obr)
             || constructIOField<vector>(e.keyword(), tok, is, obr)
             || constructIOField<sphericalTensor>(e.keyword(), tok, is, obr)
             || constructIOField<symmTensor>(e.keyword(), tok, is, obr)
             || constructIOField<tensor>(e.keyword(), tok, is, obr)
            )
            {
                // Done storing field
            }
            else
            {
                FatalErrorInFunction << "Unsupported type " << e.keyword()
                    << exit(FatalError);
            }
        }
    }
}


void Foam::mappedPatchBase::write(Ostream& os) const
{
    os.writeEntry("sampleMode", sampleModeNames_[mode_]);
    os.writeEntryIfDifferent<word>("sampleWorld", word::null, sampleWorld_);
    os.writeEntryIfDifferent<word>("sampleRegion", word::null, sampleRegion_);
    os.writeEntryIfDifferent<word>("samplePatch", word::null, samplePatch_);

    if (sampleDatabasePtr_)
    {
        os.writeEntry("sampleDatabase", Switch::name(true));
        // Write database path if differing
        os.writeEntryIfDifferent<fileName>
        (
            "sampleDatabasePath",
            fileName::null,
            sampleDatabasePtr_()
        );
    }
    coupleGroup_.write(os);

    if
    (
        offsetMode_ == UNIFORM
     && offset_ == vector::zero
     && (mode_ == NEARESTPATCHFACE || mode_ == NEARESTPATCHFACEAMI)
    )
    {
        // Collocated mode. No need to write offset data
    }
    else
    {
        os.writeEntry("offsetMode", offsetModeNames_[offsetMode_]);

        switch (offsetMode_)
        {
            case UNIFORM:
            {
                os.writeEntry("offset", offset_);
                break;
            }
            case NONUNIFORM:
            {
                offsets_.writeEntry("offsets", os);
                break;
            }
            case NORMAL:
            {
                os.writeEntry("distance", distance_);
                break;
            }
        }
    }

    if (mode_ == NEARESTPATCHFACEAMI)
    {
        if (AMIPtr_)
        {
            // Use AMI to write itself. Problem: outputs:
            // - restartUncoveredSourceFace
            // - reverseTarget (instead of flipNormals)
            AMIPtr_->write(os);
        }
        if (!surfDict_.empty())
        {
            surfDict_.writeEntry(surfDict_.dictName(), os);
        }
    }
}


// ************************************************************************* //
