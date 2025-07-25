/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022,2024 OpenCFD Ltd.
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

#include "cyclicACMIPolyPatch.H"
#include "SubField.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "primitiveMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cyclicACMIPolyPatch, 0);

    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, word);
    addToRunTimeSelectionTable(polyPatch, cyclicACMIPolyPatch, dictionary);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::cyclicACMIPolyPatch::updateAreas() const
{
    const polyMesh& mesh = boundaryMesh().mesh();

    bool updated = false;

    if (!owner())
    {
        return updated;
    }

    // Check if underlying AMI up to date
    if (!mesh.upToDatePoints(AMITime_))
    {
        // This should not happen normally since resetAMI is triggered
        // by any point motion.
        FatalErrorInFunction << "Problem : AMI is up to event:"
            << AMITime_.eventNo()
            << " mesh points are up to time " << mesh.pointsInstance()
            << " patch:" << this->name()
            << exit(FatalError);
    }

    // Check if scaling enabled (and necessary)
    if
    (
        srcScalePtr_
     && (updated || prevTimeIndex_ != mesh.time().timeIndex())
    )
    {
        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::updateAreas() :"
                << " patch:" << this->name()
                << " neighbPatch:" << this->neighbPatch().name()
                << " AMITime_:" << AMITime_.eventNo()
                << " uptodate:" << mesh.upToDatePoints(AMITime_)
                << " mesh.time().timeIndex():" << mesh.time().timeIndex()
                << " prevTimeIndex_:" << prevTimeIndex_
                << endl;
        }

        if (createAMIFaces_)
        {
            WarningInFunction
                << "Topology changes and scaling currently not supported."
                << " Patch " << this->name() << endl;
        }

        const scalar t = mesh.time().timeOutputValue();

        // Note: ideally preserve src/tgtMask before clipping to tolerance ...
        srcScaledMask_ =
            min
            (
                scalar(1) - tolerance_,
                max(tolerance_, srcScalePtr_->value(t)*srcMask_)
            );


        if (!tgtScalePtr_)
        {
            tgtScalePtr_= srcScalePtr_.clone(neighbPatch());
        }

        tgtScaledMask_ =
            min
            (
                scalar(1) - tolerance_,
                max(tolerance_, tgtScalePtr_->value(t)*tgtMask_)
            );

        if (debug)
        {
            Pout<< "cyclicACMIPolyPatch::updateAreas : scaling masks"
                << " for " << name() << " mask " << gAverage(srcScaledMask_)
                << " and " << nonOverlapPatch().name()
                << " mask " << gAverage(srcScaledMask_) << endl;
        }

        // Calculate areas from the masks
        cyclicACMIPolyPatch& cpp = const_cast<cyclicACMIPolyPatch&>(*this);
        const cyclicACMIPolyPatch& nbrCpp = neighbPatch();

        cpp.scalePatchFaceAreas(*this, srcScaledMask_, thisSf_, thisNoSf_);
        cpp.scalePatchFaceAreas(nbrCpp, tgtScaledMask_, nbrSf_, nbrNoSf_);

        prevTimeIndex_ = mesh.time().timeIndex();
        AMITime_.setUpToDate();
        updated = true;
    }

    return updated;
}


bool Foam::cyclicACMIPolyPatch::upToDate(const regIOobject& io) const
{
    // Is io up to date with
    // - underlying AMI
    // - scaling
    return io.upToDate(AMITime_);
}


void Foam::cyclicACMIPolyPatch::setUpToDate(regIOobject& io) const
{
    io.setUpToDate();
}


void Foam::cyclicACMIPolyPatch::reportCoverage
(
    const word& name,
    const scalarField& weightSum
) const
{
    label nUncovered = 0;
    label nCovered = 0;
    for (const scalar sum : weightSum)
    {
        if (sum < tolerance_)
        {
            ++nUncovered;
        }
        else if (sum > scalar(1) - tolerance_)
        {
            ++nCovered;
        }
    }
    reduce(nUncovered, sumOp<label>());
    reduce(nCovered, sumOp<label>());
    label nTotal = returnReduce(weightSum.size(), sumOp<label>());

    Info<< "ACMI: Patch " << name << " uncovered/blended/covered = "
        << nUncovered << ", " << nTotal-nUncovered-nCovered
        << ", " << nCovered << endl;
}


void Foam::cyclicACMIPolyPatch::scalePatchFaceAreas
(
    const cyclicACMIPolyPatch& acmipp,
    const scalarField& mask,                // srcMask_
    const vectorList& faceArea,             // this->faceAreas();
    const vectorList& noFaceArea            // nonOverlapPatch.faceAreas()
)
{
    // Set/scale polyPatch face areas to avoid double-accounting of face areas

    const scalar maxTol = scalar(1) - tolerance_;

    const polyPatch& nonOverlapPatch = acmipp.nonOverlapPatch();
    vectorField::subField noSf = nonOverlapPatch.faceAreas();

    DebugPout
        << "rescaling non-overlap patch areas for: "
        << nonOverlapPatch.name() << endl;

    if (mask.size() != noSf.size())
    {
        WarningInFunction
            << "Inconsistent sizes for patch: " << acmipp.name()
            << " - not manipulating patches" << nl
            << " - size: " << size() << nl
            << " - non-overlap patch size: " << noSf.size() << nl
            << " - mask size: " << mask.size() << nl
            << "This is OK for decomposition but"
            << " should be considered fatal at run-time" << endl;

        return;
    }

    {
        tmp<scalarField> scale
        (
            scalar(1)
          - min
            (
                max(mask, tolerance_),
                maxTol
            )
        );

        // Scale area
        forAll(noSf, facei)
        {
            noSf[facei] = noFaceArea[facei]*scale()[facei];
        }
        // Cache scaling
        const_cast<polyPatch&>(nonOverlapPatch).areaFraction(scale);
    }

    if (!createAMIFaces_)
    {
        // Note: for topological update (createAMIFaces_ = true)
        // AMI coupled patch face areas are updated as part of the topological
        // updates, e.g. by the calls to cyclicAMIPolyPatch's setTopology and
        // initMovePoints
        DebugPout
            << "scaling coupled patch areas for: " << acmipp.name() << endl;

        // Scale the coupled patch face areas
        vectorField::subField Sf = acmipp.faceAreas();

        tmp<scalarField> scale(max(tolerance_, mask));

        // Scale area
        forAll(Sf, facei)
        {
            Sf[facei] = faceArea[facei]*scale()[facei];
        }
        // Cache scaling
        const_cast<cyclicACMIPolyPatch&>(acmipp).areaFraction(scale);

        // Re-normalise the weights since the effect of overlap is already
        // accounted for in the area
        auto& weights = const_cast<scalarListList&>(acmipp.weights());
        auto& weightsSum = const_cast<scalarField&>(acmipp.weightsSum());
        forAll(weights, i)
        {
            scalarList& wghts = weights[i];
            if (wghts.size())
            {
                scalar& sum = weightsSum[i];

                for (scalar& w : wghts)
                {
                    w /= sum;
                }
                sum = 1.0;
            }
        }
    }

    const polyMesh& mesh = boundaryMesh().mesh();

    // Recompute the cell volumes adjacent to the patches since the cells with
    // duplicate faces are only closed after the duplicate faces have been
    // scaled.

    primitiveMeshTools::updateCellCentresAndVols
    (
        mesh,
        mesh.faceCentres(),
        mesh.faceAreas(),                       // already scaled
        uniqueSort(acmipp.faceCells()), // unique cells only
        mesh.cells(),
        const_cast<vectorField&>(mesh.cellCentres()),
        const_cast<scalarField&>(mesh.cellVolumes())
    );
}


void Foam::cyclicACMIPolyPatch::resetAMI() const
{
    resetAMI(boundaryMesh().mesh().points());
}


void Foam::cyclicACMIPolyPatch::resetAMI(const UList<point>& points) const
{
    if (!owner())
    {
        return;
    }

    const polyPatch& nonOverlapPatch = this->nonOverlapPatch();

    DebugPout
        << "cyclicACMIPolyPatch::resetAMI : recalculating weights"
        << " for " << name() << " and " << nonOverlapPatch.name()
        << endl;

    const polyMesh& mesh = boundaryMesh().mesh();

    // At this point we want face geometry but not cell geometry since we want
    // correct the face area on duplicate baffles before calculating the cell
    // centres and volumes.
    if (!mesh.hasFaceAreas())
    {
        FatalErrorInFunction
            << "primitiveMesh must already have face geometry"
            << abort(FatalError);
    }

    // Calculate the AMI using partial face-area-weighted. This leaves
    // the weights as fractions of local areas (sum(weights) = 1 means
    // face is fully covered)
    cyclicAMIPolyPatch::resetAMI(points);

    const AMIPatchToPatchInterpolation& AMI = this->AMI();

    // Output some statistics
    reportCoverage("source", AMI.srcWeightsSum());
    reportCoverage("target", AMI.tgtWeightsSum());

    // Set the mask fields
    // Note:
    // - assumes that the non-overlap patches are decomposed using the same
    //   decomposition as the coupled patches (per side)
    srcMask_ = clamp(AMI.srcWeightsSum(), zero_one{});
    tgtMask_ = clamp(AMI.tgtWeightsSum(), zero_one{});

    if (debug)
    {
        Pout<< "resetAMI" << endl;
        {
            const cyclicACMIPolyPatch& patch = *this;
            Pout<< "patch:" << patch.name() << " size:" << patch.size()
                << " non-overlap patch: " << patch.nonOverlapPatch().name()
                << " size:" << patch.nonOverlapPatch().size()
                << endl;
        }
        {
            const cyclicACMIPolyPatch& patch = this->neighbPatch();
            Pout<< "patch:" << patch.name() << " size:" << patch.size()
                << " non-overlap patch: " << patch.nonOverlapPatch().name()
                << " size:" << patch.nonOverlapPatch().size()
                << endl;
        }
    }
}


void Foam::cyclicACMIPolyPatch::scalePatchFaceAreas()
{
    if (!owner() || !canResetAMI())
    {
        return;
    }

    const polyPatch& nonOverlapPatch = this->nonOverlapPatch();
    const cyclicACMIPolyPatch& nbrPatch = this->neighbPatch();
    const polyPatch& nbrNonOverlapPatch = nbrPatch.nonOverlapPatch();

    if (srcScalePtr_)
    {
        // Use primitivePatch::faceAreas() to avoid need for additional copy?

        // Save overlap geometry for later scaling
        thisSf_ = this->faceAreas();
        thisNoSf_ = nonOverlapPatch.faceAreas();
        nbrSf_ = nbrPatch.faceAreas();
        nbrNoSf_ = nbrNonOverlapPatch.faceAreas();
    }

    // In-place scale the polyPatch face areas

    // Note
    // - using primitivePatch face areas since these are based on the raw
    //   point locations (not affected by ACMI scaling)
    scalePatchFaceAreas
    (
        *this,
        srcMask_,       // unscaled mask
        this->primitivePatch::faceAreas(),
        nonOverlapPatch.primitivePatch::faceAreas()
    );
    scalePatchFaceAreas
    (
        nbrPatch,
        tgtMask_,       // unscaled mask
        nbrPatch.primitivePatch::faceAreas(),
        nbrNonOverlapPatch.primitivePatch::faceAreas()
    );

    // Mark current AMI as up to date with points
    boundaryMesh().mesh().setUpToDatePoints(AMITime_);

    if (debug)
    {
        const auto& acmipp = *this;
        const auto& mesh = boundaryMesh().mesh();
        const auto& vols = mesh.cellVolumes();

        Info<< "cyclicACMI PP: " << acmipp.name()
            << " V:" << sum(scalarField(vols, acmipp.faceCells()))
            << nl
            << "cyclicACMI N-O PP: " << nonOverlapPatch.name()
            << " V:" << sum(scalarField(vols, nonOverlapPatch.faceCells()))
            << nl
            << "cyclicACMI NBR PP: " << nbrPatch.name()
            << " V:" << sum(scalarField(vols, nbrPatch.faceCells()))
            << nl
            << "cyclicACMI NBR N-O PP: " << nbrNonOverlapPatch.name()
            << " V:" << sum(scalarField(vols, nbrNonOverlapPatch.faceCells()))
            << nl
            << "cyclicACMI PP+N-O AREA: "
            << sum(faceAreas() + nonOverlapPatch.faceAreas()) << nl
            << "cyclicACMI NBR PP+N-O AREA: "
            << sum(nbrPatch.faceAreas() + nbrNonOverlapPatch.faceAreas())
            << endl;
    }
}


void Foam::cyclicACMIPolyPatch::initGeometry(PstreamBuffers& pBufs)
{
    DebugPout << "cyclicACMIPolyPatch::initGeometry : " << name() << endl;

    // Note: calculates transformation and triggers face centre calculation
    cyclicAMIPolyPatch::initGeometry(pBufs);

    // Initialise the AMI early to make sure we adapt the face areas before the
    // cell centre calculation gets triggered.
    if (!createAMIFaces_ && canResetAMI())
    {
        resetAMI();
    }

    scalePatchFaceAreas();
}


void Foam::cyclicACMIPolyPatch::calcGeometry(PstreamBuffers& pBufs)
{
    DebugPout << "cyclicACMIPolyPatch::calcGeometry : " << name() << endl;

    cyclicAMIPolyPatch::calcGeometry(pBufs);
}


void Foam::cyclicACMIPolyPatch::initMovePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugPout<< "cyclicACMIPolyPatch::initMovePoints : " << name() << endl;

    // Note: calculates transformation and triggers face centre calculation
    cyclicAMIPolyPatch::initMovePoints(pBufs, p);

    if (!createAMIFaces_ && canResetAMI())
    {
        resetAMI();
    }

    scalePatchFaceAreas();
}


void Foam::cyclicACMIPolyPatch::movePoints
(
    PstreamBuffers& pBufs,
    const pointField& p
)
{
    DebugPout << "cyclicACMIPolyPatch::movePoints : " << name() << endl;

    // When topology is changing, this will scale the duplicate AMI faces
    cyclicAMIPolyPatch::movePoints(pBufs, p);
}


void Foam::cyclicACMIPolyPatch::initUpdateMesh(PstreamBuffers& pBufs)
{
    DebugPout << "cyclicACMIPolyPatch::initUpdateMesh : " << name() << endl;

    cyclicAMIPolyPatch::initUpdateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::updateMesh(PstreamBuffers& pBufs)
{
    DebugPout << "cyclicACMIPolyPatch::updateMesh : " << name() << endl;

    cyclicAMIPolyPatch::updateMesh(pBufs);
}


void Foam::cyclicACMIPolyPatch::clearGeom()
{
    DebugPout << "cyclicACMIPolyPatch::clearGeom : " << name() << endl;

    cyclicAMIPolyPatch::clearGeom();
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::srcMask() const
{
    if (srcScalePtr_)
    {
        // Make sure areas are up-to-date
        updateAreas();

        return srcScaledMask_;
    }
    else
    {
        return srcMask_;
    }
}


const Foam::scalarField& Foam::cyclicACMIPolyPatch::tgtMask() const
{
    if (tgtScalePtr_)
    {
        // Make sure areas are up-to-date
        updateAreas();

        return tgtScaledMask_;
    }
    else
    {
        return tgtMask_;
    }
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * //

Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const label size,
    const label start,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const transformType transform,
    const word& defaultAMIMethod
)
:
    cyclicAMIPolyPatch
    (
        name,
        size,
        start,
        index,
        bm,
        patchType,
        transform,
        defaultAMIMethod
    ),
    nonOverlapPatchName_(),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIPtr_->setRequireMatch(false);

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const word& name,
    const dictionary& dict,
    const label index,
    const polyBoundaryMesh& bm,
    const word& patchType,
    const word& defaultAMIMethod
)
:
    cyclicAMIPolyPatch(name, dict, index, bm, patchType, defaultAMIMethod),
    nonOverlapPatchName_(dict.get<word>("nonOverlapPatch")),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    srcScalePtr_(PatchFunction1<scalar>::NewIfPresent(*this, "scale", dict)),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIPtr_->setRequireMatch(false);

    if (nonOverlapPatchName_ == name)
    {
        FatalIOErrorInFunction(dict)
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name
            << exit(FatalIOError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm
)
:
    cyclicAMIPolyPatch(pp, bm),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    srcScalePtr_(pp.srcScalePtr_.clone(*this)),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIPtr_->setRequireMatch(false);

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const label newSize,
    const label newStart,
    const word& nbrPatchName,
    const word& nonOverlapPatchName
)
:
    cyclicAMIPolyPatch(pp, bm, index, newSize, newStart, nbrPatchName),
    nonOverlapPatchName_(nonOverlapPatchName),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    srcScalePtr_(pp.srcScalePtr_.clone(*this)),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIPtr_->setRequireMatch(false);

    if (nonOverlapPatchName_ == name())
    {
        FatalErrorInFunction
            << "Non-overlapping patch name " << nonOverlapPatchName_
            << " cannot be the same as this patch " << name()
            << exit(FatalError);
    }

    // Non-overlapping patch might not be valid yet so cannot determine
    // associated patchID
}


Foam::cyclicACMIPolyPatch::cyclicACMIPolyPatch
(
    const cyclicACMIPolyPatch& pp,
    const polyBoundaryMesh& bm,
    const label index,
    const labelUList& mapAddressing,
    const label newStart
)
:
    cyclicAMIPolyPatch(pp, bm, index, mapAddressing, newStart),
    nonOverlapPatchName_(pp.nonOverlapPatchName_),
    nonOverlapPatchID_(-1),
    srcMask_(),
    tgtMask_(),
    srcScalePtr_(pp.srcScalePtr_.clone(*this)),
    AMITime_
    (
        IOobject
        (
            "AMITime",
            boundaryMesh().mesh().pointsInstance(),
            boundaryMesh().mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        dimensionedScalar("time", dimTime, -GREAT)
    ),
    prevTimeIndex_(-1)
{
    AMIPtr_->setRequireMatch(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cyclicACMIPolyPatch::newInternalProcFaces
(
    label& newFaces,
    label& newProcFaces
) const
{
    const List<labelList>& addSourceFaces = AMI().srcAddress();
    const scalarField& fMask = srcMask();

    // Add new faces as many weights for AMI
    forAll(addSourceFaces, faceI)
    {
        if (fMask[faceI] > tolerance_)
        {
            const labelList& nbrFaceIs = addSourceFaces[faceI];

            forAll(nbrFaceIs, j)
            {
                label nbrFaceI = nbrFaceIs[j];

                if (nbrFaceI < neighbPatch().size())
                {
                    // local faces
                    newFaces++;
                }
                else
                {
                    // Proc faces
                    newProcFaces++;
                }
            }
        }
    }
}


Foam::refPtr<Foam::labelListList>
Foam::cyclicACMIPolyPatch::mapCollocatedFaces() const
{
    const scalarField& fMask = srcMask();
    const labelListList& srcFaces = AMI().srcAddress();
    labelListList dOverFaces;

    dOverFaces.setSize(srcFaces.size());
    forAll(dOverFaces, faceI)
    {
        if (fMask[faceI] > tolerance_)
        {
            dOverFaces[faceI].setSize(srcFaces[faceI].size());

            forAll(dOverFaces[faceI], subFaceI)
            {
                dOverFaces[faceI][subFaceI] = srcFaces[faceI][subFaceI];
            }
        }
    }
    return refPtr<labelListList>(new labelListList(dOverFaces));
}


const Foam::cyclicACMIPolyPatch& Foam::cyclicACMIPolyPatch::neighbPatch() const
{
    const polyPatch& pp = this->boundaryMesh()[neighbPatchID()];

    // Bit of checking now we know neighbour patch
    if (!owner() && srcScalePtr_)
    {
        WarningInFunction
            << "Ignoring \"scale\" setting in slave patch " << name()
            << endl;
        srcScalePtr_.clear();
        tgtScalePtr_.clear();
    }

    return refCast<const cyclicACMIPolyPatch>(pp);
}


Foam::label Foam::cyclicACMIPolyPatch::nonOverlapPatchID() const
{
    if (nonOverlapPatchID_ == -1)
    {
        nonOverlapPatchID_ =
            this->boundaryMesh().findPatchID(nonOverlapPatchName_);

        if (nonOverlapPatchID_ == -1)
        {
            FatalErrorInFunction
                << "Illegal non-overlapping patch name " << nonOverlapPatchName_
                << nl << "Valid patch names are "
                << this->boundaryMesh().names()
                << exit(FatalError);
        }

        if (nonOverlapPatchID_ < index())
        {
            FatalErrorInFunction
                << "Boundary ordering error: " << type()
                << " patch must be defined prior to its non-overlapping patch"
                << nl
                << type() << " patch: " << name() << ", ID:" << index() << nl
                << "Non-overlap patch: " << nonOverlapPatchName_
                << ", ID:" << nonOverlapPatchID_ << nl
                << exit(FatalError);
        }

        const polyPatch& noPp = this->boundaryMesh()[nonOverlapPatchID_];

        bool ok = true;

        if (size() == noPp.size())
        {
            const scalarField magSf(mag(faceAreas()));
            const scalarField noMagSf(mag(noPp.faceAreas()));

            forAll(magSf, facei)
            {
                scalar ratio = mag(magSf[facei]/(noMagSf[facei] + ROOTVSMALL));

                if (ratio - 1 > tolerance_)
                {
                    ok = false;
                    break;
                }
            }
        }
        else
        {
            ok = false;
        }

        if (!ok)
        {
            FatalErrorInFunction
                << "Inconsistent ACMI patches " << name() << " and "
                << noPp.name() << ".  Patches should have identical topology"
                << exit(FatalError);
        }
    }

    return nonOverlapPatchID_;
}


void Foam::cyclicACMIPolyPatch::initOrder
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp
) const
{
    cyclicAMIPolyPatch::initOrder(pBufs, pp);
}


bool Foam::cyclicACMIPolyPatch::order
(
    PstreamBuffers& pBufs,
    const primitivePatch& pp,
    labelList& faceMap,
    labelList& rotation
) const
{
    return cyclicAMIPolyPatch::order(pBufs, pp, faceMap, rotation);
}


void Foam::cyclicACMIPolyPatch::write(Ostream& os) const
{
    cyclicAMIPolyPatch::write(os);

    os.writeEntry("nonOverlapPatch", nonOverlapPatchName_);

    if (owner() && srcScalePtr_)
    {
        srcScalePtr_->writeData(os);
    }
}


// ************************************************************************* //
