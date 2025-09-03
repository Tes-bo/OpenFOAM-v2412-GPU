/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020,2022,2024 OpenCFD Ltd.
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

#include "faceAreaWeightAMI2D.H"
#include "profiling.H"
#include "addToRunTimeSelectionTable.H"
#include "triangle2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceAreaWeightAMI2D, 0);
    addToRunTimeSelectionTable(AMIInterpolation, faceAreaWeightAMI2D, dict);
    addToRunTimeSelectionTable
    (
        AMIInterpolation,
        faceAreaWeightAMI2D,
        component
    );
}

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::faceAreaWeightAMI2D::writeNoMatch
(
    const label srcFacei,
    const labelList& tgtFaceCandidates,
    const boundBox& srcFaceBb
) const
{
    Info<< "NO MATCH for source face " << srcFacei << endl;
    {
        const auto& src = this->srcPatch();
        const auto& tgt = this->tgtPatch(); // might be the extended patch!

        OFstream os("no_match_" + Foam::name(srcFacei) + ".obj");

        const pointField& srcPoints = src.points();
        const pointField& tgtPoints = tgt.points();

        label np = 0;

        // Write source face
        const face& srcF = src[srcFacei];
        string faceStr = "f";
        for (const label pointi : srcF)
        {
            const point& p = srcPoints[pointi];
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
            ++np;
            faceStr += " " + Foam::name(np);
        }
        os << faceStr.c_str() << " " << (np - srcF.size() + 1) << nl;

        // Write target faces as lines
        for (const label tgtFacei : tgtFaceCandidates)
        {
            const face& tgtF = tgt[tgtFacei];
            forAll(tgtF, pointi)
            {
                const point& p = tgtPoints[tgtF[pointi]];
                os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
                ++np;
                if (pointi)
                {
                    os << "l " << np-1 << " " << np << nl;
                }
            }
            os << "l " << (np - tgtF.size() + 1) << " " << np << nl;
        }
    }

    {
        OFstream os("no_match_" + Foam::name(srcFacei) + "_bb.obj");

        const pointField points(srcFaceBb.points());
        for (const point& p : points)
        {
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << endl;
        }
        os  << "l 1 2" << nl;
        os  << "l 2 4" << nl;
        os  << "l 4 3" << nl;
        os  << "l 3 1" << nl;
        os  << "l 5 6" << nl;
        os  << "l 6 8" << nl;
        os  << "l 8 7" << nl;
        os  << "l 7 5" << nl;
        os  << "l 5 1" << nl;
        os  << "l 6 2" << nl;
        os  << "l 8 4" << nl;
        os  << "l 7 3" << nl;
    }
}


void Foam::faceAreaWeightAMI2D::storeInterArea
(
    const label srcFacei,
    const label tgtFacei,
    DynamicList<label>& srcAddr,
    DynamicList<scalar>& srcWght,
    DynamicList<vector>& srcCtr,
    DynamicList<label>& tgtAddr,
    DynamicList<scalar>& tgtWght
) const
{
    addProfiling(ami, "faceAreaWeightAMI2D::calcInterArea");

    // Quick reject if either face has zero area/too far away/wrong orientation
    if (!isCandidate(srcFacei, tgtFacei))
    {
        return;
    }

    const auto& srcPatch = this->srcPatch();
    const auto& tgtPatch = this->tgtPatch();

    const pointField& srcPoints = srcPatch.points();
    const pointField& tgtPoints = tgtPatch.points();

    const auto& srcTris = srcTris_[srcFacei];
    const auto& tgtTris = tgtTris_[tgtFacei];

    const auto& srcFaceNormals = srcPatch.faceNormals();

    //triangle2D::debug = 1;

    scalar area = 0;
    vector centroid = Zero;

    for (const auto& tris : srcTris)
    {
        const vector& origin = srcPoints[tris[0]];
        const vector p10(srcPoints[tris[1]] - origin);
        const vector p20(srcPoints[tris[2]] - origin);
        const vector axis1(p10/(mag(p10) + ROOTVSMALL));
        const vector axis2(srcFaceNormals[srcFacei]^axis1);

        triangle2D s
        (
            vector2D(0, 0),
            vector2D(axis1&p10, axis2&p10),
            vector2D(axis1&p20, axis2&p20)
        );

        for (const auto& trit : tgtTris)
        {
            // Triangle t has opposite orientation wrt triangle s
            triangle2D t
            (
                tgtPoints[trit[0]] - origin,
                tgtPoints[trit[2]] - origin,
                tgtPoints[trit[1]] - origin,
                axis1,
                axis2
            );

            scalar da = 0;
            vector2D c(Zero);

            if (t.snapClosePoints(s) == 3)
            {
                c = s.centre();
                da = mag(s.area());
            }
            else
            {
                s.interArea(t, c, da);
            }

            area += da;
            centroid += da*(origin + c.x()*axis1 + c.y()*axis2);
        }
    }

    //triangle2D::debug = 0;

    if (area/srcMagSf_[srcFacei] > triangle2D::relTol)
    {
        centroid /= area + ROOTVSMALL;

        srcAddr.append(tgtFacei);
        srcWght.append(area);
        srcCtr.append(centroid);

        tgtAddr.append(srcFacei);
        tgtWght.append(area);
    }
}


Foam::labelList Foam::faceAreaWeightAMI2D::overlappingTgtFaces
(
    const AABBTree<face>& tree,
    const List<boundBox>& tgtFaceBbs,
    const boundBox& srcFaceBb
) const
{
    labelHashSet faceIds(16);

    const auto& treeBb = tree.boundBoxes();
    const auto& treeAddr = tree.addressing();

    forAll(treeBb, boxi)
    {
        const auto& tbb = treeBb[boxi];

        if (srcFaceBb.overlaps(tbb))
        {
            const auto& boxAddr = treeAddr[boxi];

            for (const auto& tgtFacei : boxAddr)
            {
                if (srcFaceBb.overlaps(tgtFaceBbs[tgtFacei]))
                {
                    faceIds.insert(tgtFacei);
                }
            }
        }
    }

    return faceIds.toc();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceAreaWeightAMI2D::faceAreaWeightAMI2D
(
    const dictionary& dict,
    const bool reverseTarget
)
:
    advancingFrontAMI(dict, reverseTarget),
    Cbb_(dict.getCheckOrDefault<scalar>("Cbb", 0.1, scalarMinMax::ge(SMALL)))
{}


Foam::faceAreaWeightAMI2D::faceAreaWeightAMI2D
(
    const bool requireMatch,
    const bool reverseTarget,
    const scalar lowWeightCorrection,
    const faceAreaIntersect::triangulationMode triMode,
    const bool restartUncoveredSourceFace
)
:
    advancingFrontAMI
    (
        requireMatch,
        reverseTarget,
        lowWeightCorrection,
        triMode
    ),
    Cbb_(0.1)
{}


Foam::faceAreaWeightAMI2D::faceAreaWeightAMI2D
(
    const faceAreaWeightAMI2D& ami
)
:
    advancingFrontAMI(ami),
    Cbb_(ami.Cbb_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceAreaWeightAMI2D::calculate
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

    addProfiling(ami, "faceAreaWeightAMI2D::calculate");

    advancingFrontAMI::calculate(srcPatch, tgtPatch, surfPtr);

    const auto& src = this->srcPatch();
    const auto& tgt = this->tgtPatch(); // might be the extended patch!

    bool validSize = true;

    // Check that patch sizes are valid
    if (!src.size())
    {
        validSize = false;
    }
    else if (!tgt.size())
    {
        WarningInFunction
            << src.size() << " source faces but no target faces" << endl;

        validSize = false;
    }

    srcCentroids_.setSize(srcAddress_.size());

    // Temporary storage for addressing and weights
    List<DynamicList<label>> tgtAddr(tgt.size());
    List<DynamicList<scalar>> tgtWght(tgt.size());

    // Find an approximate face length scale
    const scalar lf(Cbb_*Foam::sqrt(gAverage(srcMagSf_)));

    if (validSize)
    {
        // Create the tgt tree
        const bool equalBinSize = true;
        const label maxLevel = 10;
        const label minBinSize = 4;
        AABBTree<face> tree
        (
            tgt,
            tgt.points(),
            equalBinSize,
            maxLevel,
            minBinSize
        );

        const auto& tgtPoints = tgt.points();
        List<boundBox> tgtFaceBbs(tgt.size());
        forAll(tgt, facei)
        {
            tgtFaceBbs[facei] = boundBox(tgtPoints, tgt[facei], false);
        }

        DynamicList<label> nonOverlapFaces;

        const auto& srcPoints = src.points();

        forAll(src, srcFacei)
        {
            const face& srcFace = src[srcFacei];

            treeBoundBox srcFaceBb(srcPoints, srcFace);
            srcFaceBb.grow(lf);

            const labelList tgtFaces
            (
                overlappingTgtFaces(tree, tgtFaceBbs, srcFaceBb)
            );

            DynamicList<label> srcAddr(tgtFaces.size());
            DynamicList<scalar> srcWght(tgtFaces.size());
            DynamicList<point> srcCtr(tgtFaces.size());

            for (const label tgtFacei : tgtFaces)
            {
                storeInterArea
                (
                    srcFacei,
                    tgtFacei,
                    srcAddr,
                    srcWght,
                    srcCtr,
                    tgtAddr[tgtFacei],
                    tgtWght[tgtFacei]
                );
            }

            if (mustMatchFaces() && srcAddr.empty())
            {
                if (debug) writeNoMatch(srcFacei, tgtFaces, srcFaceBb);

                // FatalErrorInFunction
                //    << "Unable to find match for source face " << srcFace
                //    << exit(FatalError);
            }

            srcAddress_[srcFacei].transfer(srcAddr);
            srcWeights_[srcFacei].transfer(srcWght);
            srcCentroids_[srcFacei].transfer(srcCtr);
        }

        srcNonOverlap_.transfer(nonOverlapFaces);

        if (debug && !srcNonOverlap_.empty())
        {
            Pout<< "    AMI: " << srcNonOverlap_.size()
                << " non-overlap faces identified"
                << endl;
        }
    }

    // Transfer data to persistent storage

    forAll(tgtAddr, i)
    {
        tgtAddress_[i].transfer(tgtAddr[i]);
        tgtWeights_[i].transfer(tgtWght[i]);
    }

    if (distributed())
    {
        const label myRank = UPstream::myProcNo(comm_);
        // Allocate unique tag for all comms
        const int oldTag = UPstream::incrMsgType();

        const primitivePatch& srcPatch0 = this->srcPatch0();
        const primitivePatch& tgtPatch0 = this->tgtPatch0();

        // Create global indexing for each original patch
        const globalIndex globalSrcFaces(srcPatch0.size(), comm_);
        const globalIndex globalTgtFaces(tgtPatch0.size(), comm_);

        for (labelList& addressing : srcAddress_)
        {
            for (label& addr : addressing)
            {
                addr = extendedTgtFaceIDs_[addr];
            }
        }

        for (labelList& addressing : tgtAddress_)
        {
            globalSrcFaces.inplaceToGlobal(myRank, addressing);
        }

        // Send data back to originating procs. Note that contributions
        // from different processors get added (ListOps::appendEqOp)

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>::null(),
            tgtPatch0.size(),
            extendedTgtMapPtr_->constructMap(),
            false,                      // has flip
            extendedTgtMapPtr_->subMap(),
            false,                      // has flip
            tgtAddress_,
            labelList(),
            ListOps::appendEqOp<label>(),
            flipOp(),                   // flip operation
            UPstream::msgType()+77431,
            comm_
        );

        mapDistributeBase::distribute
        (
            Pstream::commsTypes::nonBlocking,
            List<labelPair>::null(),
            tgtPatch0.size(),
            extendedTgtMapPtr_->constructMap(),
            false,
            extendedTgtMapPtr_->subMap(),
            false,
            tgtWeights_,
            scalarList(),
            ListOps::appendEqOp<scalar>(),
            flipOp(),                   // flip operation
            UPstream::msgType()+77432,
            comm_
        );

        // Note: using patch face areas calculated by the AMI method
        extendedTgtMapPtr_->reverseDistribute(tgtPatch0.size(), tgtMagSf_);

        // Cache maps and reset addresses
        List<Map<label>> cMapSrc;
        srcMapPtr_.reset
        (
            new mapDistribute
            (
                globalSrcFaces,
                tgtAddress_,
                cMapSrc,
                UPstream::msgType()+77433,
                comm_
            )
        );

        List<Map<label>> cMapTgt;
        tgtMapPtr_.reset
        (
            new mapDistribute
            (
                globalTgtFaces,
                srcAddress_,
                cMapTgt,
                UPstream::msgType()+77434,
                comm_
            )
        );

        // Reset tag
        UPstream::msgType(oldTag);
    }

    // Convert the weights from areas to normalised values
    normaliseWeights(requireMatch_, true);

    nonConformalCorrection();

    upToDate_ = true;

    return upToDate_;
}


void Foam::faceAreaWeightAMI2D::write(Ostream& os) const
{
    advancingFrontAMI::write(os);
    os.writeEntryIfDifferent<scalar>("Cbb", 0.1, Cbb_);
}


// ************************************************************************* //
