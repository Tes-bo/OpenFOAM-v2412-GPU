/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "nearestFaceAMI.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nearestFaceAMI, 0);
    addToRunTimeSelectionTable(AMIInterpolation, nearestFaceAMI, dict);
    addToRunTimeSelectionTable(AMIInterpolation, nearestFaceAMI, component);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::autoPtr<Foam::mapDistribute> Foam::nearestFaceAMI::calcFaceMap
(
    const List<nearestAndDist>& localInfo,
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch
) const
{
    // Generate the list of processor bounding boxes for tgtPatch
    List<boundBox> procBbs(Pstream::nProcs());
    procBbs[Pstream::myProcNo()] =
        boundBox(tgtPatch.points(), tgtPatch.meshPoints(), true);
    Pstream::allGatherList(procBbs);

    // Identify which of my local src faces intersect with each processor
    // tgtPatch bb within the current match's search distance
    const pointField& srcCcs = srcPatch.faceCentres();
    List<DynamicList<label>> dynSendMap(Pstream::nProcs());

    forAll(localInfo, srcFacei)
    {
        // Test local srcPatch face centres against remote processor tgtPatch bb
        // using distance from initial pass

        const scalar r2 = localInfo[srcFacei].second();

        forAll(procBbs, proci)
        {
            if (proci != Pstream::myProcNo())
            {
                if (procBbs[proci].overlaps(srcCcs[srcFacei], r2))
                {
                    dynSendMap[proci].append(srcFacei);
                }
            }
        }
    }

    // Convert dynamicList to labelList
    labelListList sendMap(Pstream::nProcs());
    forAll(sendMap, proci)
    {
        sendMap[proci].transfer(dynSendMap[proci]);

        if (debug)
        {
            Pout<< "send map - to proc " << proci << " sending "
                << sendMap[proci].size() << " elements" << endl;
        }
    }

    return autoPtr<mapDistribute>::New(std::move(sendMap));
}


Foam::autoPtr<Foam::mapDistribute> Foam::nearestFaceAMI::calcDistributed
(
    const primitivePatch& src,
    const primitivePatch& tgt,
    labelListList& srcToTgtAddr,
    scalarListList& srcToTgtWght
) const
{
    autoPtr<indexedOctree<treeType>> tgtTreePtr;
    if (tgt.size())
    {
        tgtTreePtr = this->createTree(tgt);
    }

    // Create global indexing for tgtPatch
    globalIndex globalTgtCells(tgt.size());


    const label myRank = UPstream::myProcNo(comm_);


    // First pass
    // ==========
    // For each srcPatch face, determine local match on tgtPatch

    // Identify local nearest matches
    pointField srcCcs(src.faceCentres());

    List<nearestAndDist> localInfo(src.size());
    if (tgt.size())
    {
        const auto& tgtTree = tgtTreePtr();

        forAll(srcCcs, srcCelli)
        {
            const point& srcCc = srcCcs[srcCelli];

            pointIndexHit& test = localInfo[srcCelli].first();
            test = tgtTree.findNearest(srcCc, GREAT);

            if (test.hit())
            {
                // With a search radius2 of GREAT all cells should receive a hit
                localInfo[srcCelli].second() = test.point().distSqr(srcCc);
                test.setIndex(globalTgtCells.toGlobal(myRank, test.index()));
            }
        }
    }
    else
    {
        // No local tgtPatch faces - initialise nearest distance for all
        // srcPatch faces to GREAT so that they [should be] set by remote
        // tgtPatch faces
        for (auto& info : localInfo)
        {
            info.second() = GREAT;
        }
    }


    // Second pass
    // ===========
    // Determine remote matches

    // Map returns labels of src patch faces to send to each proc
    autoPtr<mapDistribute> mapPtr = calcFaceMap(localInfo, src, tgt);
    mapDistribute& map = mapPtr();

    List<nearestAndDist> remoteInfo(localInfo);
    map.distribute(remoteInfo);

    // Note: re-using srcCcs
    map.distribute(srcCcs);

    if (tgt.size())
    {
        const auto& tgtTree = tgtTreePtr();

        // Test remote srcPatch faces against local tgtPatch faces
        nearestAndDist testInfo;
        pointIndexHit& test = testInfo.first();
        forAll(remoteInfo, i)
        {
            test = tgtTree.findNearest(srcCcs[i], remoteInfo[i].second());
            if (test.hit())
            {
                test.setIndex(globalTgtCells.toGlobal(myRank, test.index()));
                testInfo.second() = test.point().distSqr(srcCcs[i]);
                nearestEqOp()(remoteInfo[i], testInfo);
            }
        }
    }

    // Send back to originating processor. Choose best if sent to multiple
    // processors. Note that afterwards all unused entries have the unique
    // value nearestZero (distance < 0). This is used later on to see if
    // the sample was sent to any processor.
    const nearestAndDist nearestZero(pointIndexHit(), -GREAT);
    mapDistributeBase::distribute
    (
        Pstream::commsTypes::nonBlocking,
        List<labelPair>::null(),
        src.size(),
        map.constructMap(),
        map.constructHasFlip(),
        map.subMap(),
        map.subHasFlip(),
        remoteInfo,
        nearestZero,
        nearestEqOp(),
        identityOp(),       // No flipping
        UPstream::msgType(),
        comm_
    );


    // Third pass
    // ==========
    // Combine local and remote info and filter out any connections that are
    // further away than threshold distance squared

    srcToTgtAddr.setSize(src.size());
    srcToTgtWght.setSize(src.size());
    forAll(srcToTgtAddr, srcFacei)
    {
        nearestEqOp()(localInfo[srcFacei], remoteInfo[srcFacei]);
        if (localInfo[srcFacei].second() < maxDistance2_)
        {
            const label tgtFacei = localInfo[srcFacei].first().index();
            srcToTgtAddr[srcFacei] = labelList(1, tgtFacei);
            srcToTgtWght[srcFacei] = scalarList(1, 1.0);
        }
    }

    List<Map<label>> cMap;
    return autoPtr<mapDistribute>::New
    (
        globalTgtCells,
        srcToTgtAddr,
        cMap,
        UPstream::msgType(),
        comm_
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nearestFaceAMI::nearestFaceAMI
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    AMIInterpolation(dict, reverseTarget),
    maxDistance2_(dict.getOrDefault<scalar>("maxDistance2", GREAT))
{}


Foam::nearestFaceAMI::nearestFaceAMI
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection
)
:
    AMIInterpolation(requireMatch, reverseTarget, lowWeightCorrection),
    maxDistance2_(GREAT)
{}


Foam::nearestFaceAMI::nearestFaceAMI(const nearestFaceAMI& ami)
:
    AMIInterpolation(ami),
    maxDistance2_(ami.maxDistance2_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::nearestFaceAMI::calculate
(
    const primitivePatch& srcPatch,
    const primitivePatch& tgtPatch,
    const autoPtr<searchableSurface>& surfPtr
)
{
    if (upToDate_)
    {
        return false;
    }

    AMIInterpolation::calculate(srcPatch, tgtPatch, surfPtr);

    const auto& src = this->srcPatch0();
    const auto& tgt = this->tgtPatch0();

    // Set face area magnitudes
    srcMagSf_ = mag(src.faceAreas());
    tgtMagSf_ = mag(tgt.faceAreas());

    // TODO: implement symmetric calculation controls; assume yes for now
    bool symmetric_ = true;

    if (this->distributed())
    {
        tgtMapPtr_ =
            calcDistributed
            (
                src,
                tgt,
                srcAddress_,
                srcWeights_
            );

        if (symmetric_)
        {
            srcMapPtr_ =
                calcDistributed
                (
                    tgt,
                    src,
                    tgtAddress_,
                    tgtWeights_
                );
            }
    }
    else
    {
        srcAddress_.setSize(src.size());
        srcWeights_.setSize(src.size());

        if (symmetric_)
        {
            tgtAddress_.setSize(tgt.size());
            tgtWeights_.setSize(tgt.size());
        }

        const pointField& srcCcs = src.faceCentres();
        const pointField& tgtCcs = tgt.faceCentres();

        const auto tgtTreePtr = this->createTree(tgtPatch);
        const auto& tgtTree = tgtTreePtr();

        forAll(srcCcs, srcFacei)
        {
            const point& srcCc = srcCcs[srcFacei];
            const pointIndexHit hit = tgtTree.findNearest(srcCc, GREAT);

            if
            (
                hit.hit()
             && (magSqr(srcCc - tgtCcs[hit.index()]) < maxDistance2_)
            )
            {
                label tgtFacei = hit.index();
                srcAddress_[srcFacei] = labelList(1, tgtFacei);
                srcWeights_[srcFacei] = scalarList(1, 1.0);

                if (symmetric_)
                {
                    tgtAddress_[tgtFacei] = labelList(1, srcFacei);
                    tgtWeights_[tgtFacei] = scalarList(1, 1.0);
                }
            }
            else
            {
                if (debug)
                {
                    WarningInFunction
                        << "Unable to find target face for source face "
                        << srcFacei << endl;
                }
            }
        }

        if (symmetric_)
        {
            const auto srcTreePtr = this->createTree(srcPatch);
            const auto& srcTree = srcTreePtr();

            // Check that all source cells have connections and populate any
            // missing entries
            forAll(tgtWeights_, tgtCelli)
            {
                if (tgtAddress_[tgtCelli].empty())
                {
                    const point& tgtCc = tgtCcs[tgtCelli];
                    pointIndexHit hit = srcTree.findNearest(tgtCc, GREAT);

                    if
                    (
                        hit.hit()
                     && (magSqr(tgtCc - srcCcs[hit.index()]) < maxDistance2_)
                    )
                    {
                        tgtAddress_[tgtCelli] = labelList(1, hit.index());
                        tgtWeights_[tgtCelli] = scalarList(1, 1.0);
                    }
                }
                else
                {
                    if (debug)
                    {
                        WarningInFunction
                            << "Unable to find source face for target face "
                            << tgtCelli << endl;
                    }
                }
            }
        }
    }

    srcWeightsSum_.setSize(srcWeights_.size(), 1);
    tgtWeightsSum_.setSize(tgtWeights_.size(), 1);

    upToDate_ = true;

    return upToDate_;
}


void Foam::nearestFaceAMI::write(Ostream& os) const
{
    AMIInterpolation::write(os);
    os.writeEntryIfDifferent<scalar>("maxDistance2", GREAT, maxDistance2_);
}


// ************************************************************************* //
