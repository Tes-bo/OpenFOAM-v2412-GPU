/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "turbulentDFSEMInletFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "momentOfInertia.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::label Foam::turbulentDFSEMInletFvPatchVectorField::seedIterMax_ = 1000;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::turbulentDFSEMInletFvPatchVectorField::writeEddyOBJ() const
{
    {
        // Output the bounding box
        OFstream os(db().time().path()/"eddyBox.obj");

        const polyPatch& pp = this->patch().patch();
        const labelList& boundaryPoints = pp.boundaryPoints();
        const pointField& localPoints = pp.localPoints();

        const vector offset(patchNormal_*maxSigmaX_);
        forAll(boundaryPoints, i)
        {
            point p = localPoints[boundaryPoints[i]];
            p += offset;
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
        }

        forAll(boundaryPoints, i)
        {
            point p = localPoints[boundaryPoints[i]];
            p -= offset;
            os  << "v " << p.x() << " " << p.y() << " " << p.z() << nl;
        }
    }

    {
        const Time& time = db().time();
        OFstream os
        (
            time.path()/"eddies_" + Foam::name(time.timeIndex()) + ".obj"
        );

        label pointOffset = 0;
        forAll(eddies_, eddyI)
        {
            const eddy& e = eddies_[eddyI];
            pointOffset += e.writeSurfaceOBJ(pointOffset, patchNormal_, os);
        }
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::writeLumleyCoeffs() const
{
    // Output list of xi vs eta

    OFstream os(db().time().path()/"lumley_interpolated.out");

    os  << "# xi" << token::TAB << "eta" << endl;

    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    forAll(R, faceI)
    {
        // Normalised anisotropy tensor
        const symmTensor devR(dev(R[faceI]/(tr(R[faceI]))));

        // Second tensor invariant
        const scalar ii = min(0, invariantII(devR));

        // Third tensor invariant
        const scalar iii = invariantIII(devR);

        // xi, eta
        // See Pope - characterization of Reynolds-stress anisotropy
        const scalar xi = cbrt(0.5*iii);
        const scalar eta = sqrt(-ii/3.0);
        os  << xi << token::TAB << eta << token::TAB
            << ii << token::TAB << iii << endl;
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialisePatch()
{
    const vectorField nf(patch().nf());

    // Patch normal points into domain
    patchNormal_ = -gAverage(nf);

    // Check that patch is planar
    const scalar error = max(magSqr(patchNormal_ + nf));

    if (error > SMALL)
    {
        WarningInFunction
            << "Patch " << patch().name() << " is not planar"
            << endl;
    }

    patchNormal_ /= mag(patchNormal_) + ROOTVSMALL;


    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    // Triangulate the patch faces and create addressing
    {
        label nTris = 0;
        for (const face& f : patch)
        {
            nTris += f.nTriangles();
        }

        DynamicList<labelledTri> dynTriFace(nTris);
        DynamicList<face> tris(8);  // work array

        forAll(patch, facei)
        {
            const face& f = patch[facei];

            tris.clear();
            f.triangles(points, tris);

            for (const auto& t : tris)
            {
                dynTriFace.emplace_back(t[0], t[1], t[2], facei);
            }
        }

        // Transfer to persistent storage
        triFace_.transfer(dynTriFace);
    }


    const label myProci = UPstream::myProcNo();
    const label numProc = UPstream::nProcs();

    // Calculate the cumulative triangle weights
    {
        triCumulativeMagSf_.resize_nocopy(triFace_.size()+1);

        auto iter = triCumulativeMagSf_.begin();

        // Set zero value at the start of the tri area/weight list
        scalar patchArea = 0;
        *iter++ = patchArea;

        // Calculate cumulative and total area
        for (const auto& t : triFace_)
        {
            patchArea += t.mag(points);
            *iter++ = patchArea;
        }

        sumTriMagSf_.resize_nocopy(numProc+1);
        sumTriMagSf_[0] = 0;

        {
            scalarList::subList slice(sumTriMagSf_, numProc, 1);
            slice[myProci] = patchArea;
            Pstream::allGatherList(slice);
        }

        // Convert to cumulative sum of areas per proc
        for (label i = 1; i < sumTriMagSf_.size(); ++i)
        {
            sumTriMagSf_[i] += sumTriMagSf_[i-1];
        }
    }

    // Global patch area (over all processors)
    patchArea_ = sumTriMagSf_.back();

    // Local patch bounds (this processor)
    patchBounds_ = boundBox(patch.localPoints(), false);
    patchBounds_.inflate(0.1);

    // Determine if all eddies spawned from a single processor
    singleProc_ = returnReduceOr
    (
        patch.size() == returnReduce(patch.size(), sumOp<label>())
    );
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddyBox()
{
    const scalarField& magSf = patch().magSf();

    const scalarField L(L_->value(db().time().timeOutputValue())/Lref_);

    // (PCF:Eq. 14)
    const scalarField cellDx(max(Foam::sqrt(magSf), 2/patch().deltaCoeffs()));

    // Inialise eddy box extents
    forAll(*this, faceI)
    {
        scalar& s = sigmax_[faceI];

        // Average length scale (SST:Eq. 24)
        // Personal communication regarding (PCR:Eq. 14)
        //  - the min operator in Eq. 14 is a typo, and should be a max operator
        s = min(mag(L[faceI]), kappa_*delta_);
        s = max(s, nCellPerEddy_*cellDx[faceI]);
    }

    // Maximum extent across all processors
    maxSigmaX_ = gMax(sigmax_);

    // Eddy box volume
    v0_ = 2*gSum(magSf)*maxSigmaX_;

    if (debug)
    {
        Info<< "Patch: " << patch().patch().name() << " eddy box:" << nl
            << "    volume    : " << v0_ << nl
            << "    maxSigmaX : " << maxSigmaX_ << nl
            << endl;
    }
}


Foam::pointIndexHit Foam::turbulentDFSEMInletFvPatchVectorField::setNewPosition
(
    const bool global
)
{
    const polyPatch& patch = this->patch().patch();
    const pointField& points = patch.points();

    label triI = -1;

    if (global)
    {
        const scalar areaFraction =
            rndGen_.globalPosition<scalar>(0, patchArea_);

        // Determine which processor to use
        label proci = 0;
        forAllReverse(sumTriMagSf_, i)
        {
            if (areaFraction >= sumTriMagSf_[i])
            {
                proci = i;
                break;
            }
        }

        if (UPstream::myProcNo() == proci)
        {
            // Find corresponding decomposed face triangle
            triI = 0;
            const scalar offset = sumTriMagSf_[proci];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    triI = i;
                    break;
                }
            }
        }
    }
    else
    {
        // Find corresponding decomposed face triangle on local processor
        triI = 0;
        const scalar maxAreaLimit = triCumulativeMagSf_.back();
        const scalar areaFraction = rndGen_.position<scalar>(0, maxAreaLimit);

        forAllReverse(triCumulativeMagSf_, i)
        {
            if (areaFraction > triCumulativeMagSf_[i])
            {
                triI = i;
                break;
            }
        }
    }


    if (triI >= 0)
    {
        return pointIndexHit
        (
            true,
            // Find random point in triangle
            triFace_[triI].tri(points).randomPoint(rndGen_),
            triFace_[triI].index()
        );
    }

    // No hit
    return pointIndexHit(false, vector::max, -1);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::initialiseEddies()
{
    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    DynamicList<eddy> eddies(size());

    // Initialise eddy properties
    scalar sumVolEddy = 0;
    scalar sumVolEddyAllProc = 0;

    while (sumVolEddyAllProc/v0_ < d_)
    {
        bool search = true;
        label iter = 0;

        while (search && iter++ < seedIterMax_)
        {
            // Get new parallel consistent position
            pointIndexHit pos(setNewPosition(true));
            const label patchFaceI = pos.index();

            // Note: only 1 processor will pick up this face
            if (patchFaceI != -1)
            {
                eddy e
                (
                    patchFaceI,
                    pos.hitPoint(),
                    rndGen_.position<scalar>(-maxSigmaX_, maxSigmaX_),
                    sigmax_[patchFaceI],
                    R[patchFaceI],
                    rndGen_
                );

                // If eddy valid, patchFaceI is non-zero
                if (e.patchFaceI() != -1)
                {
                    eddies.append(e);
                    sumVolEddy += e.volume();
                    search = false;
                }
            }
            // else eddy on remote processor

            reduce(search, andOp<bool>());
        }


        sumVolEddyAllProc = returnReduce(sumVolEddy, sumOp<scalar>());
    }
    eddies_.transfer(eddies);

    nEddy_ = eddies_.size();

    if (debug)
    {
        Pout<< "Patch:" << patch().patch().name();

        if (Pstream::parRun())
        {
            Pout<< " processor:" << Pstream::myProcNo();
        }

        Pout<< " seeded:" << nEddy_ << " eddies" << endl;
    }

    reduce(nEddy_, sumOp<label>());

    if (nEddy_ > 0)
    {
        Info<< "Turbulent DFSEM patch: " << patch().name()
            << " seeded " << nEddy_ << " eddies with total volume "
            << sumVolEddyAllProc
            << endl;
    }
    else
    {
        WarningInFunction
            << "Patch: " << patch().patch().name()
            << " on field " << internalField().name()
            << ": No eddies seeded - please check your set-up"
            << endl;
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::convectEddies
(
    const vector& UBulk,
    const scalar deltaT
)
{
    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    // Note: all operations applied to local processor only

    label nRecycled = 0;

    forAll(eddies_, eddyI)
    {
        eddy& e = eddies_[eddyI];
        e.move(deltaT*(UBulk & patchNormal_));

        const scalar position0 = e.x();

        // Check to see if eddy has exited downstream box plane
        if (position0 > maxSigmaX_)
        {
            bool search = true;
            label iter = 0;

            while (search && iter++ < seedIterMax_)
            {
                // Spawn new eddy with new random properties (intensity etc)
                pointIndexHit pos(setNewPosition(false));
                const label patchFaceI = pos.index();

                e = eddy
                    (
                        patchFaceI,
                        pos.hitPoint(),
                        position0 - floor(position0/maxSigmaX_)*maxSigmaX_,
                        sigmax_[patchFaceI],
                        R[patchFaceI],
                        rndGen_
                    );

                if (e.patchFaceI() != -1)
                {
                    search = false;
                }
            }

            nRecycled++;
        }
    }

    if (debug)
    {
        reduce(nRecycled, sumOp<label>());

        if (nRecycled)
        {
            Info<< "Patch: " << patch().patch().name()
                << " recycled " << nRecycled << " eddies"
                << endl;
        }
    }
}


Foam::vector Foam::turbulentDFSEMInletFvPatchVectorField::uPrimeEddy
(
    const List<eddy>& eddies,
    const point& patchFaceCf
) const
{
    vector uPrime(Zero);

    forAll(eddies, k)
    {
        const eddy& e = eddies[k];
        uPrime += e.uPrime(patchFaceCf, patchNormal_);
    }

    return uPrime;
}


void Foam::turbulentDFSEMInletFvPatchVectorField::calcOverlappingProcEddies
(
    List<List<eddy>>& overlappingEddies
) const
{
    const int oldTag = UPstream::incrMsgType();

    List<boundBox> patchBBs(Pstream::nProcs());
    patchBBs[Pstream::myProcNo()] = patchBounds_;
    Pstream::allGatherList(patchBBs);

    // Per processor indices into all segments to send
    List<DynamicList<label>> sendMap(UPstream::nProcs());

    // Collect overlapping eddies
    forAll(eddies_, i)
    {
        const eddy& e = eddies_[i];

        // Eddy bounds
        const point x(e.position(patchNormal_));
        boundBox ebb(e.bounds());
        ebb.min() += x;
        ebb.max() += x;

        forAll(patchBBs, proci)
        {
            // Not including intersection with local patch
            if (proci != Pstream::myProcNo())
            {
                if (ebb.overlaps(patchBBs[proci]))
                {
                    sendMap[proci].push_back(i);
                }
            }
        }
    }


    PstreamBuffers pBufs;

    forAll(sendMap, proci)
    {
        if (proci != Pstream::myProcNo() && !sendMap[proci].empty())
        {
            UOPstream os(proci, pBufs);

            os << UIndirectList<eddy>(eddies_, sendMap[proci]);
        }
    }

    pBufs.finishedSends();

    for (const int proci : pBufs.allProcs())
    {
        if (proci != Pstream::myProcNo() && pBufs.recvDataCount(proci))
        {
            UIPstream is(proci, pBufs);

            is >> overlappingEddies[proci];
        }
    }

    UPstream::msgType(oldTag);  // Restore tag
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    U_(nullptr),
    R_(nullptr),
    L_(nullptr),
    delta_(1.0),
    d_(1.0),
    kappa_(0.41),
    Uref_(1.0),
    Lref_(1.0),
    scale_(1.0),
    m_(0.5),
    nCellPerEddy_(5),

    patchArea_(-1),
    triFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    patchNormal_(Zero),
    patchBounds_(),

    eddies_(Zero),
    v0_(Zero),
    rndGen_(Pstream::myProcNo()),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    nEddy_(0),
    curTimeIndex_(-1),
    singleProc_(false),
    writeEddies_(false)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    U_(ptf.U_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    L_(ptf.L_.clone(patch().patch())),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    m_(ptf.m_),
    nCellPerEddy_(ptf.nCellPerEddy_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    patchNormal_(ptf.patchNormal_),
    patchBounds_(ptf.patchBounds_),

    eddies_(ptf.eddies_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_, mapper),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(ptf.nEddy_),
    curTimeIndex_(-1),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    U_(PatchFunction1<vector>::New(patch().patch(), "U", dict)),
    R_(PatchFunction1<symmTensor>::New(patch().patch(), "R", dict)),
    L_(PatchFunction1<scalar>::New(patch().patch(), "L", dict)),
    delta_(dict.getCheck<scalar>("delta", scalarMinMax::ge(0))),
    d_(dict.getCheckOrDefault<scalar>("d", 1, scalarMinMax::ge(SMALL))),
    kappa_(dict.getCheckOrDefault<scalar>("kappa", 0.41, scalarMinMax::ge(0))),
    Uref_(dict.getCheckOrDefault<scalar>("Uref", 1, scalarMinMax::ge(SMALL))),
    Lref_(dict.getCheckOrDefault<scalar>("Lref", 1, scalarMinMax::ge(SMALL))),
    scale_(dict.getCheckOrDefault<scalar>("scale", 1, scalarMinMax::ge(0))),
    m_(dict.getCheckOrDefault<scalar>("m", 0.5, scalarMinMax::ge(0))),
    nCellPerEddy_(dict.getOrDefault<label>("nCellPerEddy", 5)),

    patchArea_(-1),
    triFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero),
    patchNormal_(Zero),
    patchBounds_(),

    eddies_(),
    v0_(Zero),
    rndGen_(),
    sigmax_(size(), Zero),
    maxSigmaX_(Zero),
    nEddy_(0),
    curTimeIndex_(-1),
    singleProc_(false),
    writeEddies_(dict.getOrDefault("writeEddies", false))
{
    eddy::debug = debug;

    const scalar t = db().time().timeOutputValue();
    const symmTensorField R(R_->value(t)/sqr(Uref_));

    checkStresses(R);
}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    U_(ptf.U_.clone(patch().patch())),
    R_(ptf.R_.clone(patch().patch())),
    L_(ptf.L_.clone(patch().patch())),
    delta_(ptf.delta_),
    d_(ptf.d_),
    kappa_(ptf.kappa_),
    Uref_(ptf.Uref_),
    Lref_(ptf.Lref_),
    scale_(ptf.scale_),
    m_(ptf.m_),
    nCellPerEddy_(ptf.nCellPerEddy_),

    patchArea_(ptf.patchArea_),
    triFace_(ptf.triFace_),
    triCumulativeMagSf_(ptf.triCumulativeMagSf_),
    sumTriMagSf_(ptf.sumTriMagSf_),
    patchNormal_(ptf.patchNormal_),
    patchBounds_(ptf.patchBounds_),

    eddies_(ptf.eddies_),
    v0_(ptf.v0_),
    rndGen_(ptf.rndGen_),
    sigmax_(ptf.sigmax_),
    maxSigmaX_(ptf.maxSigmaX_),
    nEddy_(ptf.nEddy_),
    curTimeIndex_(-1),
    singleProc_(ptf.singleProc_),
    writeEddies_(ptf.writeEddies_)
{}


Foam::turbulentDFSEMInletFvPatchVectorField::
turbulentDFSEMInletFvPatchVectorField
(
    const turbulentDFSEMInletFvPatchVectorField& ptf
)
:
    turbulentDFSEMInletFvPatchVectorField(ptf, ptf.internalField())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentDFSEMInletFvPatchVectorField::checkStresses
(
    const symmTensorField& R
)
{
    constexpr label maxDiffs = 5;
    label nDiffs = 0;

    // (S:Eq. 4a-4c)
    forAll(R, i)
    {
        bool diff = false;

        if (maxDiffs < nDiffs)
        {
            Info<< "More than " << maxDiffs << " times"
                << " Reynolds-stress realizability checks failed."
                << " Skipping further comparisons." << endl;
            return;
        }

        const symmTensor& r = R[i];

        if (r.xx() < 0)
        {
            WarningInFunction
                << "Reynolds stress " << r << " at index " << i
                << " does not obey the constraint: Rxx >= 0"
                << endl;
            diff = true;
        }

        if ((r.xx()*r.yy() - sqr(r.xy())) < 0)
        {
            WarningInFunction
                << "Reynolds stress " << r << " at index " << i
                << " does not obey the constraint: Rxx*Ryy - sqr(Rxy) >= 0"
                << endl;
            diff = true;
        }

        if (det(r) < 0)
        {
            WarningInFunction
                << "Reynolds stress " << r << " at index " << i
                << " does not obey the constraint: det(R) >= 0"
                << endl;
            diff = true;
        }

        if (diff)
        {
            ++nDiffs;
        }
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::checkStresses
(
    const scalarField& R
)
{
    if (min(R) <= 0)
    {
        FatalErrorInFunction
            << "Reynolds stresses contain at least one "
            << "nonpositive element. min(R) = " << min(R)
            << exit(FatalError);
    }
}


void Foam::turbulentDFSEMInletFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<vector>::autoMap(m);

    if (U_)
    {
        U_->autoMap(m);
    }
    if (R_)
    {
        R_->autoMap(m);
    }
    if (L_)
    {
        L_->autoMap(m);
    }

    sigmax_.autoMap(m);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<vector>::rmap(ptf, addr);

    const auto& dfsemptf =
        refCast<const turbulentDFSEMInletFvPatchVectorField>(ptf);

    if (U_)
    {
        U_->rmap(dfsemptf.U_(), addr);
    }
    if (R_)
    {
        R_->rmap(dfsemptf.R_(), addr);
    }
    if (L_)
    {
        L_->rmap(dfsemptf.L_(), addr);
    }

    sigmax_.rmap(dfsemptf.sigmax_, addr);
}


void Foam::turbulentDFSEMInletFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ == -1)
    {
        initialisePatch();

        initialiseEddyBox();

        initialiseEddies();
    }


    if (curTimeIndex_ != db().time().timeIndex())
    {
        tmp<vectorField> UMean =
            U_->value(db().time().timeOutputValue())/Uref_;

        // (PCR:p. 522)
        const vector UBulk
        (
            gSum(UMean()*patch().magSf())
           /(gSum(patch().magSf()) + ROOTVSMALL)
        );

        // Move eddies using bulk velocity
        const scalar deltaT = db().time().deltaTValue();
        convectEddies(UBulk, deltaT);

        // Set mean velocity
        vectorField& U = *this;
        U = UMean;

        // Apply second part of normalisation coefficient
        const scalar c =
            scale_*Foam::pow(10*v0_, m_)/Foam::sqrt(scalar(nEddy_));

        // In parallel, need to collect all eddies that will interact with
        // local faces

        const pointField& Cf = patch().Cf();

        if (singleProc_ || !Pstream::parRun())
        {
            forAll(U, faceI)
            {
                U[faceI] += c*uPrimeEddy(eddies_, Cf[faceI]);
            }
        }
        else
        {
            // Process local eddy contributions
            forAll(U, faceI)
            {
                U[faceI] += c*uPrimeEddy(eddies_, Cf[faceI]);
            }

            // Add contributions from overlapping eddies
            List<List<eddy>> overlappingEddies(Pstream::nProcs());
            calcOverlappingProcEddies(overlappingEddies);

            for (const List<eddy>& eddies : overlappingEddies)
            {
                if (eddies.size())
                {
                    forAll(U, faceI)
                    {
                        U[faceI] += c*uPrimeEddy(eddies, Cf[faceI]);
                    }
                }
            }
        }

        // Re-scale to ensure correct flow rate
        const scalar fCorr =
            gSum((UBulk & patchNormal_)*patch().magSf())
           /gSum(U & -patch().Sf());

        U *= fCorr;

        curTimeIndex_ = db().time().timeIndex();

        if (writeEddies_)
        {
            writeEddyOBJ();
        }

        if (debug)
        {
            Info<< "Magnitude of bulk velocity: " << UBulk << endl;

            Info<< "Number of eddies: "
                << returnReduce(eddies_.size(), sumOp<label>())
                << endl;

            Info<< "Patch:" << patch().patch().name()
                << " min/max(U):" << gMin(U) << ", " << gMax(U)
                << endl;

            if (db().time().writeTime())
            {
                writeLumleyCoeffs();
            }
        }
    }

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::turbulentDFSEMInletFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("delta", delta_);
    os.writeEntryIfDifferent<scalar>("d", 1.0, d_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("Uref", 1.0, Uref_);
    os.writeEntryIfDifferent<scalar>("Lref", 1.0, Lref_);
    os.writeEntryIfDifferent<scalar>("scale", 1.0, scale_);
    os.writeEntryIfDifferent<scalar>("m", 0.5, m_);
    os.writeEntryIfDifferent<label>("nCellPerEddy", 5, nCellPerEddy_);
    os.writeEntryIfDifferent("writeEddies", false, writeEddies_);
    if (U_)
    {
        U_->writeData(os);
    }
    if (R_)
    {
        R_->writeData(os);
    }
    if (L_)
    {
        L_->writeData(os);
    }
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       turbulentDFSEMInletFvPatchVectorField
   );
}


// ************************************************************************* //
