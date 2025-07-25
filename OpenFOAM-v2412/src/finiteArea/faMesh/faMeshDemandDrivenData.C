/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
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

#include "faMesh.H"
#include "faMeshLduAddressing.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "fac.H"
#include "cartesianCS.H"
#include "scalarMatrices.H"
#include "processorFaPatch.H"
#include "processorFaPatchFields.H"
#include "emptyFaPatchFields.H"
#include "wedgeFaPatch.H"
#include "triangle.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Define an area-weighted normal for three points (defining a triangle)
// (p0, p1, p2) are the base, first, second points respectively
//
// From the original Tukovic code:
//
//   vector n = (d1^d2)/mag(d1^d2);
//   scalar sinAlpha = mag(d1^d2)/(mag(d1)*mag(d2));
//   scalar w = sinAlpha/(mag(d1)*mag(d2));
//
// ie, normal weighted by area, sine angle and inverse distance squared.
// - area : larger weight for larger areas
// - sin  : lower weight for narrow angles (eg, shards)
// - inv distance squared : lower weights for distant points
//
// The above refactored, with 0.5 for area:
//
//   (d1 ^ d2) / (2 * magSqr(d1) * magSqr(d2))

static inline vector areaInvDistSqrWeightedNormal
(
    const vector& a,
    const vector& b
)
{
    // TBD: Use volatile to avoid aggressive branch optimization
    const scalar s(2*magSqr(a)*magSqr(b));

    return s < ROOTVSMALL ? Zero : (a ^ b) / s;
}


// The area normal for the face dual (around base-point)
// described by the right/left edge points and the centre point
//
// The adjustment for 1/2 edge point (the dual point) is done internally
static inline vector areaInvDistSqrWeightedNormalDualEdge
(
    const point& basePoint,
    const point& rightPoint,
    const point& leftPoint,
    const point& centrePoint
)
{
    const vector mid(centrePoint - basePoint);

    return
    (
        areaInvDistSqrWeightedNormal
        (
            0.5*(rightPoint - basePoint),   // vector to right 1/2 edge
            mid
        )
      + areaInvDistSqrWeightedNormal
        (
            mid,
            0.5*(leftPoint - basePoint)     // vector to left 1/2 edge
        )
    );
}


// Weighted area normal for the face triangle (around base-point)
// to the edge points and the centre point.
// Used for example for area-weighted edge normals
// ---------------------------------------------------------------------------
//         own          |
//          *           | Don't bother trying to split the triangle into
//         /.\          | sub-triangles. Planar anyhow and skew weighting
//        / . \         | is less relevant here.
//       /  .  \        |
//      /   .   \       | Note: need negative for neighbour side orientation.
//  e0 *----:----* e1   |
// ---------------------------------------------------------------------------

static inline vector areaInvDistSqrWeightedNormalFaceTriangle
(
    const linePointRef& edgeLine,
    const point& faceCentre
)
{
    return
    (
        areaInvDistSqrWeightedNormal
        (
            (edgeLine.a() - faceCentre),  // From centre to right edge
            (edgeLine.b() - faceCentre)   // From centre to left edge
        )
    );
}


// Simple area normal calculation for specified edge vector and face centre
// The face centre comes last since it is less accurate
static inline vector areaNormalFaceTriangle
(
    const linePointRef& edgeLine,
    const point& faceCentre
)
{
    return triPointRef::areaNormal(edgeLine.a(), edgeLine.b(), faceCentre);
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Calculate transform tensor with reference vector (unitAxis1)
// and direction vector (axis2).
//
// This is nearly identical to the meshTools axesRotation
// with an E3_E1 transformation with the following exceptions:
//
// - axis1 (e3 == unitAxis1): is already normalized (unit vector)
// - axis2 (e1 == dirn): no difference
// - transformation is row-vectors, not column-vectors
static inline tensor rotation_e3e1
(
    const vector& unitAxis1,
    vector dirn
)
{
    dirn.removeCollinear(unitAxis1);
    dirn.normalise();

    // Set row vectors
    return tensor
    (
        dirn,
        (unitAxis1^dirn),
        unitAxis1
    );
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// A bitSet (size patch nPoints()) with boundary points marked
static Foam::bitSet markupBoundaryPoints(const uindirectPrimitivePatch& p)
{
    // Initially all unmarked
    bitSet markPoints(p.nPoints());
    for (const edge& e : p.boundaryEdges())
    {
        // Mark boundary points
        markPoints.set(e.first());
        markPoints.set(e.second());
    }

    return markPoints;
}


} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::faMesh::calcLduAddressing() const
{
    DebugInFunction
        << "Calculating addressing" << endl;

    if (lduPtr_)
    {
        FatalErrorInFunction
            << "lduPtr_ already allocated"
            << abort(FatalError);
    }

    lduPtr_ = std::make_unique<faMeshLduAddressing>(*this);
}


void Foam::faMesh::calcPatchStarts() const
{
    DebugInFunction
        << "Calculating patch starts" << endl;

    if (patchStartsPtr_)
    {
        FatalErrorInFunction
            << "patchStarts already allocated"
            << abort(FatalError);
    }

    patchStartsPtr_ = std::make_unique<labelList>(boundary().patchStarts());
}


void Foam::faMesh::calcWhichPatchFaces() const
{
    // Usually need both together
    if (polyPatchFacesPtr_ || polyPatchIdsPtr_)
    {
        FatalErrorInFunction
            << "Already allocated polyPatchFaces/polyPatchIds"
            << abort(FatalError);
    }

    const polyBoundaryMesh& pbm = mesh().boundaryMesh();

    polyPatchFacesPtr_ = std::make_unique<List<labelPair>>
    (
        pbm.whichPatchFace(faceLabels_)
    );

    labelHashSet ids;

    // Extract patch ids from (patch, facei) tuples
    for (const labelPair& tup : *polyPatchFacesPtr_)
    {
        ids.insert(tup.first());
    }

    ids.erase(-1);  // Without internal faces (patchi == -1)

    // parSync
    Foam::reduce
    (
        ids,
        bitOrOp<labelHashSet>(),
        UPstream::msgType(),
        this->comm()
    );

    polyPatchIdsPtr_ = std::make_unique<labelList>(ids.sortedToc());
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

//- Calculate the 'Le' vector from faceCentre to edge centre
//  using the edge normal to correct for curvature
//
//  Normalise and rescaled to the edge length
static inline vector calcLeVector
(
    const point& faceCentre,
    const linePointRef& edgeLine,
    const vector& edgeNormal   // (unit or area normal)
)
{
    const vector centreToEdge(edgeLine.centre() - faceCentre);

    vector leVector(edgeLine.vec() ^ edgeNormal);

    scalar s(mag(leVector));

    if (s < ROOTVSMALL)
    {
        // The calculated edgeNormal somehow degenerate and thus a
        // bad cross-product?
        // Revert to basic centre -> edge

        leVector = centreToEdge;
        leVector.removeCollinear(edgeLine.unitVec());
        s = mag(leVector);

        if (s < ROOTVSMALL)
        {
            // Unlikely that this should happen
            return Zero;
        }

        leVector *= edgeLine.mag()/s;
    }
    else
    {
        // The additional orientation is probably unnecessary
        leVector *= edgeLine.mag()/s * sign(leVector & centreToEdge);
    }

    return leVector;
}

} // End namespace Foam


namespace Foam
{

// Impose a minimum (edge) vector length
// - disallow any mag(vec) < SMALL
static inline void imposeMinVectorLength(vectorField& fld)
{
    // If it is small, no orientation information reasonably possible
    // so use a vector(1,1,1) * sqrt(1/3)*min-length

    // sqrt(1/3) = 0.5773502691896257, but slightly rounded down
    const vector minVector(vector::uniform(0.57735*SMALL));
    const scalar minLenSqr(SMALL*SMALL);

    for (vector& v : fld)
    {
        if (v.magSqr() < minLenSqr)
        {
            v = minVector;
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::vectorField> Foam::faMesh::calcRawEdgeNormals(int order) const
{
    // Return edge normals with flat boundary addressing
    auto tedgeNormals = tmp<vectorField>::New(nEdges_);
    auto& edgeNormals = tedgeNormals.ref();

    // Need face centres
    const areaVectorField& fCentres = areaCentres();

    // Also need local points
    const pointField& localPoints = points();


    if (order < 0)
    {
        // geometryOrder (-1): no communication

        // Simple (primitive) edge normal calculations.
        // These are primarly designed to avoid any communication
        // but are thus necessarily inconsistent across processor boundaries!

        WarningInFunction
            << "Using geometryOrder < 0 : "
               "simplified edge area-normals, without processor connectivity"
            << endl;

        // Internal (edge normals) - contributions from owner/neighbour
        for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
        {
            const linePointRef edgeLine(edges_[edgei].line(localPoints));

            edgeNormals[edgei] =
            (
                areaNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeOwner()[edgei]]
                )
                // NB: reversed sign (edge orientation flipped for neighbour)
              - areaNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeNeighbour()[edgei]]
                )
            );
        }

        // Boundary (edge normals). Like above, but only has owner
        for (label edgei = nInternalEdges_; edgei < nEdges_; ++edgei)
        {
            const linePointRef edgeLine(edges_[edgei].line(localPoints));

            edgeNormals[edgei] =
            (
                areaNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeOwner()[edgei]]
                )
            );
        }
    }
    else
    {
        // geometryOrder (0) - no communication
        // otherwise with communication

        if (order < 1)
        {
            WarningInFunction
                << "Using geometryOrder == 0 : "
                   "weighted edge normals, without processor connectivity"
                << endl;
        }

        // Internal (edge normals)
        // - area-weighted contributions from owner/neighbour
        for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
        {
            const linePointRef edgeLine(edges_[edgei].line(localPoints));

            edgeNormals[edgei] =
            (
                areaInvDistSqrWeightedNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeOwner()[edgei]]
                )
                // NB: reversed sign (edge orientation flipped for neighbour)
              - areaInvDistSqrWeightedNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeNeighbour()[edgei]]
                )
            );
        }

        // Boundary (edge normals). Like above, but only has owner
        for (label edgei = nInternalEdges_; edgei < nEdges_; ++edgei)
        {
            const linePointRef edgeLine(edges_[edgei].line(localPoints));

            edgeNormals[edgei] =
            (
                areaInvDistSqrWeightedNormalFaceTriangle
                (
                    edgeLine,
                    fCentres[edgeOwner()[edgei]]
                )
            );
        }
    }


    // Boundary edge corrections
    bitSet nbrBoundaryAdjust;

    // Processor-processor first (for convenience)
    if (order >= 1)
    {
        DynamicList<UPstream::Request> requests
        (
            2*(boundary().size() - boundary().nNonProcessor())
        );

        PtrList<vectorField> nbrEdgeNormals(boundary().size());

        // Prefer our own slicing into the raw list (instead of patchSlice).
        // Although patchSlice does properly handle the start() offset in
        // the presence of emptyPaPatch, we are slightly paranoia,
        // and want to avoid kicking off the patchStarts() data...

        label patchOffset = nInternalEdges_;

        forAll(boundary(), patchi)
        {
            const faPatch& fap = boundary()[patchi];

            auto edgeNorms = edgeNormals.slice(patchOffset, fap.nEdges());
            patchOffset += edgeNorms.size();

            const auto* fapp = isA<processorFaPatch>(fap);

            if (fapp)
            {
                const label nbrProci = fapp->neighbProcNo();

                if (UPstream::parRun() && !edgeNorms.empty())
                {
                    nbrEdgeNormals.emplace(patchi, edgeNorms.size());

                    // Recv accumulated weighted edge normals
                    UIPstream::read
                    (
                        requests.emplace_back(),
                        nbrProci,
                        nbrEdgeNormals[patchi]
                    );

                    // Send accumulated weighted edge normals
                    UOPstream::write
                    (
                        requests.emplace_back(),
                        nbrProci,
                        edgeNorms
                    );
                }
            }
            else if (isA<wedgeFaPatch>(fap))
            {
                // Correct wedge edges
                const auto& wedgePatch = refCast<const wedgeFaPatch>(fap);

                const vector wedgeNorm
                (
                    transform
                    (
                        wedgePatch.edgeT(),
                        wedgePatch.centreNormal()
                    ).normalise()
                );

                for (vector& e : edgeNorms)
                {
                    e.removeCollinear(wedgeNorm);
                }
            }
            // TBD:
            /// else if (isA<emptyFaPatch>(fap))
            /// {
            ///     // Ignore this boundary
            /// }
            else if (correctPatchPointNormals(patchi) && !fap.coupled())
            {
                // Neighbour correction requested
                nbrBoundaryAdjust.set(patchi);
            }
        }

        // Wait for outstanding requests
        UPstream::waitRequests(requests);

        if (!requests.empty())
        {
            // Add in received values

            patchOffset = nInternalEdges_;

            forAll(boundary(), patchi)
            {
                const faPatch& fap = boundary()[patchi];

                auto edgeNorms = edgeNormals.slice(patchOffset, fap.nEdges());
                patchOffset += edgeNorms.size();

                if (nbrEdgeNormals.set(patchi))
                {
                    const vectorField& nbrNorms = nbrEdgeNormals[patchi];

                    forAll(edgeNorms, patchEdgei)
                    {
                        edgeNorms[patchEdgei] += nbrNorms[patchEdgei];
                    }
                }
            }
        }
    }


    // Apply boundary edge corrections
    if
    (
        order >= 1
     && hasPatchPointNormalsCorrection()
     && returnReduceOr(nbrBoundaryAdjust.any())
    )
    {
        DebugInFunction
            << "Apply " << nbrBoundaryAdjust.count()
            << " boundary neighbour corrections" << nl;

        // Collect face normals, per boundary edge

        (void)this->haloFaceNormals();

        // Using our own slicing (instead of patchSlice) - see above
        label patchOffset = nInternalEdges_;

        for (const label patchi : nbrBoundaryAdjust)
        {
            const faPatch& fap = boundary()[patchi];

            auto edgeNorms = edgeNormals.slice(patchOffset, fap.nEdges());
            patchOffset += edgeNorms.size();

            if (fap.ngbPolyPatchIndex() < 0)
            {
                FatalErrorInFunction
                    << "Neighbour polyPatch index is not defined "
                    << "for faPatch " << fap.name()
                    << abort(FatalError);
            }

            // Apply the correction

            // Note from Zeljko Tukovic:
            //
            // This posibility is used for free-surface tracking
            // calculations to enforce 90 deg contact angle between
            // finite-area mesh and neighbouring polyPatch. It is very
            // important for curvature calculation to have correct normal
            // at contact line points.

            vectorField nbrNorms(this->haloFaceNormals(patchi));

            const label nPatchEdges = min(edgeNorms.size(), nbrNorms.size());

            for (label patchEdgei = 0; patchEdgei < nPatchEdges; ++patchEdgei)
            {
                edgeNorms[patchEdgei].removeCollinear(nbrNorms[patchEdgei]);
            }
        }
    }


    // Remove collinear components and normalise

    forAll(edgeNormals, edgei)
    {
        const linePointRef edgeLine(edges_[edgei].line(localPoints));

        edgeNormals[edgei].removeCollinear(edgeLine.unitVec());
        edgeNormals[edgei].normalise();
    }

    imposeMinVectorLength(edgeNormals);

    return tedgeNormals;
}


void Foam::faMesh::calcLe() const
{
    DebugInFunction
        << "Calculating local edges" << endl;

    if (LePtr_)
    {
        FatalErrorInFunction
            << "LePtr_ already allocated"
            << abort(FatalError);
    }

    LePtr_ = std::make_unique<edgeVectorField>
    (
        IOobject
        (
            "Le",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength
        // -> calculatedType()
    );
    edgeVectorField& Le = *LePtr_;

    // Need face centres
    const areaVectorField& fCentres = areaCentres();

    // Also need local points
    const pointField& localPoints = points();


    // The edgeAreaNormals _may_ use communication (depends on geometryOrder)

    {
        const edgeVectorField& edgeNormals = edgeAreaNormals();

        // Internal (edge vector)
        {
            vectorField& fld = Le.primitiveFieldRef();
            for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
            {
                fld[edgei] = calcLeVector
                (
                    fCentres[edgeOwner()[edgei]],
                    edges_[edgei].line(localPoints),
                    edgeNormals[edgei]
                );
            }
        }

        // Boundary (edge vector)
        forAll(boundary(), patchi)
        {
            const faPatch& fap = boundary()[patchi];
            vectorField& pfld = Le.boundaryFieldRef()[patchi];

            const vectorField& bndEdgeNormals =
                edgeNormals.boundaryField()[patchi];

            label edgei = fap.start();

            forAll(pfld, patchEdgei)
            {
                pfld[patchEdgei] = calcLeVector
                (
                    fCentres[edgeOwner()[edgei]],
                    edges_[edgei].line(localPoints),
                    bndEdgeNormals[patchEdgei]
                );

                ++edgei;
            }
        }
    }


    // Impose a minimum (edge) vector length

    // Internal (edge vector)
    {
        imposeMinVectorLength(Le.primitiveFieldRef());
    }

    // Boundary (edge vector)
    for (vectorField& pfld : Le.boundaryFieldRef())
    {
        imposeMinVectorLength(pfld);
    }
}


void Foam::faMesh::calcMagLe() const
{
    DebugInFunction
        << "Calculating local edge magnitudes" << endl;

    if (magLePtr_)
    {
        FatalErrorInFunction
            << "magLe() already allocated"
            << abort(FatalError);
    }

    magLePtr_ = std::make_unique<edgeScalarField>
    (
        IOobject
        (
            "magLe",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength
    );
    edgeScalarField& magLe = *magLePtr_;

    const pointField& localPoints = points();

    // Internal (edge length)
    {
        auto iter = magLe.primitiveFieldRef().begin();

        for (const edge& e : internalEdges())
        {
            *iter = e.mag(localPoints);

            // Do not allow any mag(val) < SMALL
            if (mag(*iter) < SMALL)
            {
                *iter = SMALL;
            }

            ++iter;
        }
    }

    // Boundary (edge length)
    {
        auto& bfld = magLe.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.mag(localPoints);

                // Do not allow any mag(val) < SMALL
                if (mag(*iter) < SMALL)
                {
                    *iter = SMALL;
                }

                ++iter;
            }
        }
    }
}


void Foam::faMesh::calcFaceCentres() const
{
    DebugInFunction
        << "Calculating face centres" << endl;

    if (faceCentresPtr_)
    {
        FatalErrorInFunction
            << "areaCentres already allocated"
            << abort(FatalError);
    }

    faceCentresPtr_ = std::make_unique<areaVectorField>
    (
        IOobject
        (
            "centres",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength
        // -> calculatedType()
    );
    areaVectorField& centres = *faceCentresPtr_;

    // Need local points
    const pointField& localPoints = points();


    // Internal (face centres)
    {
        if (mesh().hasFaceCentres())
        {
            // The volume mesh has faceCentres, can reuse them

            centres.primitiveFieldRef()
                = UIndirectList<vector>(mesh().faceCentres(), faceLabels());
        }
        else
        {
            // Calculate manually
            auto iter = centres.primitiveFieldRef().begin();

            for (const face& f : faces())
            {
                *iter = f.centre(localPoints);
                ++iter;
            }
        }
    }

    // Boundary (edge centres or neighbour face centres)
    {
        auto& bfld = centres.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.centre(localPoints);
                ++iter;
            }
        }
    }

    // Parallel consistency, exchange on processor patches
    if (UPstream::parRun())
    {
        centres.boundaryFieldRef().evaluateCoupled<processorFaPatch>();
    }
}


void Foam::faMesh::calcEdgeCentres() const
{
    DebugInFunction
        << "Calculating edge centres" << endl;

    if (edgeCentresPtr_)
    {
        FatalErrorInFunction
            << "edgeCentres already allocated"
            << abort(FatalError);
    }

    edgeCentresPtr_ = std::make_unique<edgeVectorField>
    (
        IOobject
        (
            "edgeCentres",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimLength
        // -> calculatedType()
    );
    edgeVectorField& centres = *edgeCentresPtr_;

    // Need local points
    const pointField& localPoints = points();


    // Internal (edge centres)
    {
        auto iter = centres.primitiveFieldRef().begin();

        for (const edge& e : internalEdges())
        {
            *iter = e.centre(localPoints);
            ++iter;
        }
    }

    // Boundary (edge centres)
    {
        auto& bfld = centres.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            auto iter = bfld[patchi].begin();

            for (const edge& e : boundary()[patchi].patchSlice(edges_))
            {
                *iter = e.centre(localPoints);
                ++iter;
            }
        }
    }
}


void Foam::faMesh::calcS() const
{
    DebugInFunction
        << "Calculating areas" << endl;

    if (SPtr_)
    {
        FatalErrorInFunction
            << "S() already allocated"
            << abort(FatalError);
    }

    SPtr_ = std::make_unique<DimensionedField<scalar, areaMesh>>
    (
        IOobject
        (
            "S",
            time().timeName(),
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimArea
    );
    auto& areas = *SPtr_;


    // No access to fvMesh::magSf(), only polyMesh::faceAreas()
    if (mesh().hasFaceAreas())
    {
        // The volume mesh has faceAreas, can reuse them
        UIndirectList<vector> meshFaceAreas(mesh().faceAreas(), faceLabels());

        auto& fld = areas.field();

        forAll(fld, facei)
        {
            fld[facei] = Foam::mag(meshFaceAreas[facei]);

            // Do not allow any mag(val) < SMALL
            if (mag(fld[facei]) < SMALL)
            {
                fld[facei] = SMALL;
            }
        }
    }
    else
    {
        // Calculate manually

        const pointField& localPoints = points();

        auto iter = areas.field().begin();

        for (const face& f : faces())
        {
            *iter = f.mag(localPoints);

            // Do not allow any mag(val) < SMALL
            if (mag(*iter) < SMALL)
            {
                *iter = SMALL;
            }

            ++iter;
        }
    }
}


void Foam::faMesh::calcFaceAreaNormals() const
{
    DebugInFunction
         << "Calculating face area normals" << endl;

    if (faceAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "faceAreaNormals already allocated"
            << abort(FatalError);
    }

    faceAreaNormalsPtr_ = std::make_unique<areaVectorField>
    (
        IOobject
        (
            "faceAreaNormals",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless
    );
    areaVectorField& faceNormals = *faceAreaNormalsPtr_;

    const pointField& localPoints = points();

    // Internal
    {
        auto& fld = faceNormals.primitiveFieldRef();

        if (mesh().hasFaceAreas())
        {
            // The volume mesh has faceAreas, can reuse them
            fld = UIndirectList<vector>(mesh().faceAreas(), faceLabels());
        }
        else
        {
            // Calculate manually

            auto iter = fld.begin();

            for (const face& f : faces())
            {
                *iter = f.areaNormal(localPoints);
                ++iter;
            }
        }

        // Make unit normals
        fld.normalise();

        imposeMinVectorLength(fld);
    }


    // Boundary - copy from edges
    {
        const auto& edgeNormalsBoundary = edgeAreaNormals().boundaryField();

        forAll(boundary(), patchi)
        {
            faceNormals.boundaryFieldRef()[patchi]
                = edgeNormalsBoundary[patchi];
        }
    }

    // Parallel consistency, exchange on processor patches
    if (UPstream::parRun())
    {
        faceNormals.boundaryFieldRef().evaluateCoupled<processorFaPatch>();
    }
}


void Foam::faMesh::calcEdgeAreaNormals() const
{
    DebugInFunction
        << "Calculating edge area normals" << endl;

    if (edgeAreaNormalsPtr_)
    {
        FatalErrorInFunction
            << "edgeAreaNormalsPtr_ already allocated"
            << abort(FatalError);
    }

    edgeAreaNormalsPtr_ = std::make_unique<edgeVectorField>
    (
        IOobject
        (
            "edgeAreaNormals",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless
        // -> calculatedType()
    );
    edgeVectorField& edgeAreaNormals = *edgeAreaNormalsPtr_;


    if (faMesh::geometryOrder() < 2)
    {
        // The edge normals with flat boundary addressing
        // (which _may_ use communication)
        vectorField edgeNormals
        (
            calcRawEdgeNormals(faMesh::geometryOrder())
        );

        // Copy internal field
        edgeAreaNormals.primitiveFieldRef()
            = vectorField::subList(edgeNormals, nInternalEdges_);

        // Using our own slicing (instead of patchSlice) - see above

        auto& bfld = edgeAreaNormals.boundaryFieldRef();

        label patchOffset = nInternalEdges_;

        forAll(boundary(), patchi)
        {
            const faPatch& fap = boundary()[patchi];

            bfld[patchi] = vectorField::subList
            (
                edgeNormals,
                bfld[patchi].size(),
                patchOffset
            );

            patchOffset += fap.nEdges();
        }

        return;
    }


    // ------------------------------------------------------------------------

    // This is the original approach using an average of the
    // point normals.  May be removed in the future (2022-09)

    // Starting from point area normals
    const vectorField& pointNormals = pointAreaNormals();

    // Also need local points
    const pointField& localPoints = points();

    // Internal edges
    {
        auto& fld = edgeAreaNormals.primitiveFieldRef();

        forAll(fld, edgei)
        {
            const edge& e = edges_[edgei];
            const linePointRef edgeLine(e.line(localPoints));

            // Average of both ends
            fld[edgei] = (pointNormals[e.first()] + pointNormals[e.second()]);

            fld[edgei].removeCollinear(edgeLine.unitVec());
            fld[edgei].normalise();
        }

        imposeMinVectorLength(fld);
    }

    // Boundary
    {
        auto& bfld = edgeAreaNormals.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            const faPatch& fap = boundary()[patchi];

            auto& pfld = bfld[patchi];

            const edgeList::subList bndEdges = fap.patchSlice(edges_);

            forAll(bndEdges, patchEdgei)
            {
                const edge& e = bndEdges[patchEdgei];
                const linePointRef edgeLine(e.line(localPoints));

                // Average of both ends
                pfld[patchEdgei] =
                    (pointNormals[e.first()] + pointNormals[e.second()]);

                pfld[patchEdgei].removeCollinear(edgeLine.unitVec());
                pfld[patchEdgei].normalise();
            }

            imposeMinVectorLength(pfld);
        }
    }
}


void Foam::faMesh::calcFaceCurvatures() const
{
    DebugInFunction
        << "Calculating face curvatures" << endl;

    if (faceCurvaturesPtr_)
    {
        FatalErrorInFunction
            << "faceCurvaturesPtr_ already allocated"
            << abort(FatalError);
    }

    faceCurvaturesPtr_ = std::make_unique<areaScalarField>
    (
        IOobject
        (
            "faceCurvatures",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless/dimLength
    );
    areaScalarField& faceCurvatures = *faceCurvaturesPtr_;


//     faceCurvatures =
//         fac::edgeIntegrate(Le()*edgeLengthCorrection())
//         &faceAreaNormals();

    areaVectorField kN(fac::edgeIntegrate(Le()*edgeLengthCorrection()));

    faceCurvatures = sign(kN&faceAreaNormals())*mag(kN);
}


void Foam::faMesh::calcEdgeTransformTensors() const
{
    DebugInFunction
        << "Calculating edge transformation tensors" << endl;

    if (edgeTransformTensorsPtr_)
    {
        FatalErrorInFunction
            << "edgeTransformTensorsPtr_ already allocated"
            << abort(FatalError);
    }

    edgeTransformTensorsPtr_ =
        std::make_unique<FieldField<Field, tensor>>(nEdges_);
    auto& edgeTransformTensors = *edgeTransformTensorsPtr_;

    // Initialize all transformation tensors
    for (label edgei = 0; edgei < nEdges_; ++edgei)
    {
        edgeTransformTensors.set(edgei, new Field<tensor>(3, tensor::I));
    }

    const areaVectorField& Nf = faceAreaNormals();
    const areaVectorField& Cf = areaCentres();

    const edgeVectorField& Ne = edgeAreaNormals();
    const edgeVectorField& Ce = edgeCentres();

    // Internal edges transformation tensors
    for (label edgei = 0; edgei < nInternalEdges_; ++edgei)
    {
        const label ownFacei = owner()[edgei];
        const label neiFacei = neighbour()[edgei];
        auto& tensors = edgeTransformTensors[edgei];

        vector edgeCtr = Ce.internalField()[edgei];

        if (skew())
        {
            edgeCtr -= skewCorrectionVectors().internalField()[edgei];
        }

        // Edge transformation tensor
        tensors[0] =
            rotation_e3e1
            (
                Ne.internalField()[edgei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Owner transformation tensor
        tensors[1] =
            rotation_e3e1
            (
                Nf.internalField()[ownFacei],
                (edgeCtr - Cf.internalField()[ownFacei])
            );

        // Neighbour transformation tensor
        tensors[2] =
            rotation_e3e1
            (
                Nf.internalField()[neiFacei],
                (Cf.internalField()[neiFacei] - edgeCtr)
            );
    }

    // Boundary edges transformation tensors
    forAll(boundary(), patchi)
    {
        const labelUList& edgeFaces = boundary()[patchi].edgeFaces();

        if (boundary()[patchi].coupled())
        {
            vectorField nbrCf(Cf.boundaryField()[patchi].patchNeighbourField());
            vectorField nbrNf(Nf.boundaryField()[patchi].patchNeighbourField());

            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];

                if (skew())
                {
                    edgeCtr -= skewCorrectionVectors()
                        .boundaryField()[patchi][bndEdgei];
                }

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] =
                    rotation_e3e1
                    (
                        nbrNf[bndEdgei],
                        (nbrCf[bndEdgei] - edgeCtr)
                    );
            }
        }
        else
        {
            forAll(edgeFaces, bndEdgei)
            {
                const label ownFacei = edgeFaces[bndEdgei];
                const label meshEdgei(boundary()[patchi].start() + bndEdgei);

                auto& tensors = edgeTransformTensors[meshEdgei];

                vector edgeCtr = Ce.boundaryField()[patchi][bndEdgei];

                if (skew())
                {
                    edgeCtr -= skewCorrectionVectors()
                        .boundaryField()[patchi][bndEdgei];
                }

                // Edge transformation tensor
                tensors[0] =
                    rotation_e3e1
                    (
                        Ne.boundaryField()[patchi][bndEdgei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Owner transformation tensor
                tensors[1] =
                    rotation_e3e1
                    (
                        Nf.internalField()[ownFacei],
                        (edgeCtr - Cf.internalField()[ownFacei])
                    );

                // Neighbour transformation tensor
                tensors[2] = tensor::I;
            }
        }
    }
}


Foam::labelList Foam::faMesh::internalPoints() const
{
    DebugInFunction
        << "Calculating internal points" << endl;

    bitSet markPoints(markupBoundaryPoints(this->patch()));
    markPoints.flip();

    return markPoints.sortedToc();
}


Foam::labelList Foam::faMesh::boundaryPoints() const
{
    DebugInFunction
        << "Calculating boundary points" << endl;

    bitSet markPoints(markupBoundaryPoints(this->patch()));

    return markPoints.sortedToc();
}


// ~~~~~~~~~~~~~~~~~~~~~~~~~
// Point normal calculations
// ~~~~~~~~~~~~~~~~~~~~~~~~~

// Revised method (general)
// ------------------------
//
// - For each patch face obtain a centre point (mathematical avg)
//   and use that to define the face dual as a pair of triangles:
//     - tri1: base-point / mid-side of right edge / face centre
//     - tri2: base-point / face centre / mid-side of left edge
//
// - Walk all face points, inserting directly into the corresponding
//   locations. No distinction between internal or boundary points (yet).
//
// Revised method (boundary correction)
// ------------------------------------
//
// - correct wedge directly, use processor patch information to exchange
//   the current summed values. [No change from original].
//
// - explicit correction of other boundaries.
//   Use the new boundary halo information for the face normals.
//   Calculate the equivalent face-point normals locally and apply
//   correction as before.

void Foam::faMesh::calcPointAreaNormals(vectorField& result) const
{
    DebugInFunction
        << "Calculating pointAreaNormals : face dual method" << endl;

    result.resize_nocopy(nPoints());
    result = Zero;

    const pointField& points = patch().localPoints();
    const faceList& faces = patch().localFaces();


    // Loop over all faces

    for (const face& f : faces)
    {
        const label nVerts(f.size());

        point centrePoint(Zero);
        for (const label fp : f)
        {
            centrePoint += points[fp];
        }
        centrePoint /= nVerts;

        for (label i = 0; i < nVerts; ++i)
        {
            const label pt0 = f.thisLabel(i);    // base

            result[pt0] +=
                areaInvDistSqrWeightedNormalDualEdge
                (
                    points[pt0],             // base
                    points[f.nextLabel(i)],  // right
                    points[f.prevLabel(i)],  // left
                    centrePoint
                );
        }
    }


    // Boundary edge corrections
    bitSet nbrBoundaryAdjust;

    forAll(boundary(), patchi)
    {
        const faPatch& fap = boundary()[patchi];

        if (isA<wedgeFaPatch>(fap))
        {
            // Correct wedge points

            const auto& wedgePatch = refCast<const wedgeFaPatch>(fap);

            const labelList& patchPoints = wedgePatch.pointLabels();

            const vector N
            (
                transform
                (
                    wedgePatch.edgeT(),
                    wedgePatch.centreNormal()
                ).normalise()
            );

            for (const label pti : patchPoints)
            {
                result[pti].removeCollinear(N);
            }

            // Axis point correction
            if (wedgePatch.axisPoint() > -1)
            {
                result[wedgePatch.axisPoint()] =
                    wedgePatch.axis()
                   *(
                       wedgePatch.axis()
                      &result[wedgePatch.axisPoint()]
                    );
            }
        }
        else if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            // Correct processor patch points
            const auto& procPatch = refCast<const processorFaPatch>(fap);

            const labelList& patchPoints = procPatch.pointLabels();
            const labelList& nbrPatchPoints = procPatch.neighbPoints();

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            // Send my values

            vectorField patchPointNormals
            (
                UIndirectList<vector>(result, patchPoints)
            );

            {
                UOPstream::write
                (
                    UPstream::commsTypes::buffered,
                    procPatch.neighbProcNo(),
                    patchPointNormals
                );
            }

            // Receive neighbour values
            patchPointNormals.resize_nocopy(nbrPatchPoints.size());

            {
                UIPstream::read
                (
                    UPstream::commsTypes::buffered,
                    procPatch.neighbProcNo(),
                    patchPointNormals
                );
            }

            for (const label pti : nonGlobalPatchPoints)
            {
                result[patchPoints[pti]] +=
                    patchPointNormals[nbrPatchPoints[pti]];
            }
        }
        // TBD:
        /// else if (isA<emptyFaPatch>(fap))
        /// {
        ///     // Ignore this boundary
        /// }
        else if (correctPatchPointNormals(patchi) && !fap.coupled())
        {
            // Neighbour correction requested
            nbrBoundaryAdjust.set(patchi);
        }
    }


    // Correct global points
    if (globalData().nGlobalPoints())
    {
        const labelList& spLabels = globalData().sharedPointLabels();
        const labelList& addr = globalData().sharedPointAddr();

        vectorField spNormals
        (
            UIndirectList<vector>(result, spLabels)
        );

        vectorField gpNormals
        (
            globalData().nGlobalPoints(),
            Zero
        );

        forAll(addr, i)
        {
            gpNormals[addr[i]] += spNormals[i];
        }

        Pstream::combineReduce(gpNormals, plusEqOp<vectorField>());

        // Extract local data
        forAll(addr, i)
        {
            spNormals[i] = gpNormals[addr[i]];
        }

        forAll(spNormals, pointI)
        {
            result[spLabels[pointI]] = spNormals[pointI];
        }
    }


    if
    (
        hasPatchPointNormalsCorrection()
     && returnReduceOr(nbrBoundaryAdjust.any())
    )
    {
        DebugInFunction
            << "Apply " << nbrBoundaryAdjust.count()
            << " boundary neighbour corrections" << nl;

        // Apply boundary points correction

        // Collect face normals as point normals

        const auto& haloNormals = this->haloFaceNormals();

        Map<vector> fpNormals(4*nBoundaryEdges());

        for (const label patchi : nbrBoundaryAdjust)
        {
            const faPatch& fap = boundary()[patchi];

            if (fap.ngbPolyPatchIndex() < 0)
            {
                FatalErrorInFunction
                    << "Neighbour polyPatch index is not defined "
                    << "for faPatch " << fap.name()
                    << abort(FatalError);
            }

            // NB: haloFaceNormals uses primitivePatch edge indexing
            for (const label edgei : fap.edgeLabels())
            {
                const edge& e = patch().edges()[edgei];

                // Halo face unitNormal at boundary edge (starts as 0)
                const vector& fnorm = haloNormals[edgei - nInternalEdges_];

                fpNormals(e.first()) += fnorm;
                fpNormals(e.second()) += fnorm;
            }
        }

        // Apply the correction

        // Note from Zeljko Tukovic:
        //
        // This posibility is used for free-surface tracking
        // calculations to enforce 90 deg contact angle between
        // finite-area mesh and neighbouring polyPatch. It is very
        // important for curvature calculation to have correct normal
        // at contact line points.

        forAllConstIters(fpNormals, iter)
        {
            const label pointi = iter.key();
            result[pointi].removeCollinear(normalised(iter.val()));
        }
    }

    result.normalise();
}


void Foam::faMesh::calcPointAreaNormalsByQuadricsFit(vectorField& result) const
{
    const labelList intPoints(internalPoints());
    const labelList bndPoints(boundaryPoints());

    const pointField& points = patch().localPoints();
    const faceList& faces = patch().localFaces();
    const labelListList& pointFaces = patch().pointFaces();

    forAll(intPoints, pointI)
    {
        label curPoint = intPoints[pointI];

        const labelHashSet faceSet(pointFaces[curPoint]);
        const labelList curFaces(faceSet.toc());

        labelHashSet pointSet;

        for (const label facei : curFaces)
        {
            const labelList& facePoints = faces[facei];
            pointSet.insert(facePoints);
        }
        pointSet.erase(curPoint);
        labelList curPoints(pointSet.toc());

        if (pointSet.size() < 5)
        {
            DebugInfo
                << "WARNING: Extending point set for fitting." << endl;

            labelHashSet faceSet(pointFaces[curPoint]);
            labelList curFaces(faceSet.toc());
            for (const label facei : curFaces)
            {
                const labelList& curFaceFaces = patch().faceFaces()[facei];
                faceSet.insert(curFaceFaces);
            }
            curFaces = faceSet.toc();

            pointSet.clear();

            for (const label facei : curFaces)
            {
                const labelList& facePoints = faces[facei];
                pointSet.insert(facePoints);
            }
            pointSet.erase(curPoint);
            curPoints = pointSet.toc();
        }

        pointField allPoints(curPoints.size());
        scalarField W(curPoints.size(), 1.0);
        for (label i=0; i<curPoints.size(); ++i)
        {
            allPoints[i] = points[curPoints[i]];
            W[i] = 1.0/magSqr(allPoints[i] - points[curPoint]);
        }

        // Transform points
        coordSystem::cartesian cs
        (
            points[curPoint],   // origin
            result[curPoint],   // axis [e3] (normalized by constructor)
            allPoints[0] - points[curPoint] // direction [e1]
        );

        for (point& p : allPoints)
        {
            p = cs.localPosition(p);
        }

        scalarRectangularMatrix M
        (
            allPoints.size(),
            5,
            0.0
        );

        for (label i = 0; i < allPoints.size(); ++i)
        {
            M[i][0] = sqr(allPoints[i].x());
            M[i][1] = sqr(allPoints[i].y());
            M[i][2] = allPoints[i].x()*allPoints[i].y();
            M[i][3] = allPoints[i].x();
            M[i][4] = allPoints[i].y();
        }

        scalarSquareMatrix MtM(5, Zero);

        for (label i = 0; i < MtM.n(); ++i)
        {
            for (label j = 0; j < MtM.m(); ++j)
            {
                for (label k = 0; k < M.n(); ++k)
                {
                    MtM[i][j] += M[k][i]*M[k][j]*W[k];
                }
            }
        }

        scalarField MtR(5, Zero);

        for (label i=0; i<MtR.size(); ++i)
        {
            for (label j=0; j<M.n(); ++j)
            {
                MtR[i] += M[j][i]*allPoints[j].z()*W[j];
            }
        }

        LUsolve(MtM, MtR);

        vector curNormal = vector(MtR[3], MtR[4], -1);

        curNormal = cs.globalVector(curNormal);

        curNormal *= sign(curNormal&result[curPoint]);

        result[curPoint] = curNormal;
    }


    forAll(boundary(), patchI)
    {
        const faPatch& fap = boundary()[patchI];

        if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(boundary()[patchI]);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList toNgbProcLsPointStarts(patchPointLabels.size(), Zero);
            vectorField toNgbProcLsPoints
            (
                10*patchPointLabels.size(),
                Zero
            );
            label nPoints = 0;

            for (label pointI=0; pointI<patchPointLabels.size(); ++pointI)
            {
                label curPoint = patchPointLabels[pointI];

                toNgbProcLsPointStarts[pointI] = nPoints;

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;

                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                for (label i=0; i<curPoints.size(); ++i)
                {
                    toNgbProcLsPoints[nPoints++] = points[curPoints[i]];
                }
            }

            toNgbProcLsPoints.setSize(nPoints);

            {
                OPstream toNeighbProc
                (
                    Pstream::commsTypes::buffered,
                    procPatch.neighbProcNo(),
                    toNgbProcLsPoints.size_bytes()
                  + toNgbProcLsPointStarts.size_bytes()
                  + 10*sizeof(label)
                );

                toNeighbProc
                    << toNgbProcLsPoints
                    << toNgbProcLsPointStarts;
            }
        }
    }

    for (const faPatch& fap : boundary())
    {
        if (Pstream::parRun() && isA<processorFaPatch>(fap))
        {
            const processorFaPatch& procPatch =
                refCast<const processorFaPatch>(fap);

            const labelList& patchPointLabels = procPatch.pointLabels();

            labelList fromNgbProcLsPointStarts(patchPointLabels.size(), Zero);
            vectorField fromNgbProcLsPoints;

            {
                IPstream fromNeighbProc
                (
                    Pstream::commsTypes::buffered,
                    procPatch.neighbProcNo(),
                    10*patchPointLabels.size()*sizeof(vector)
                  + fromNgbProcLsPointStarts.size_bytes()
                  + 10*sizeof(label)
                );

                fromNeighbProc
                    >> fromNgbProcLsPoints
                    >> fromNgbProcLsPointStarts;
            }

            const labelList& nonGlobalPatchPoints =
                procPatch.nonGlobalPatchPoints();

            forAll(nonGlobalPatchPoints, pointI)
            {
                label curPoint =
                    patchPointLabels[nonGlobalPatchPoints[pointI]];
                label curNgbPoint =
                    procPatch.neighbPoints()[nonGlobalPatchPoints[pointI]];

                labelHashSet faceSet;
                faceSet.insert(pointFaces[curPoint]);

                labelList curFaces = faceSet.toc();

                labelHashSet pointSet;

                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                label nAllPoints = curPoints.size();

                if (curNgbPoint == fromNgbProcLsPointStarts.size() - 1)
                {
                    nAllPoints +=
                        fromNgbProcLsPoints.size()
                      - fromNgbProcLsPointStarts[curNgbPoint];
                }
                else
                {
                    nAllPoints +=
                        fromNgbProcLsPointStarts[curNgbPoint + 1]
                      - fromNgbProcLsPointStarts[curNgbPoint];
                }

                vectorField allPointsExt(nAllPoints);
                label counter = 0;
                for (label i=0; i<curPoints.size(); ++i)
                {
                    allPointsExt[counter++] = points[curPoints[i]];
                }

                if (curNgbPoint == fromNgbProcLsPointStarts.size() - 1)
                {
                    for
                    (
                        label i=fromNgbProcLsPointStarts[curNgbPoint];
                        i<fromNgbProcLsPoints.size();
                        ++i
                    )
                    {
                        allPointsExt[counter++] = fromNgbProcLsPoints[i];
                    }
                }
                else
                {
                    for
                    (
                        label i=fromNgbProcLsPointStarts[curNgbPoint];
                        i<fromNgbProcLsPointStarts[curNgbPoint+1];
                        ++i
                    )
                    {
                        allPointsExt[counter++] = fromNgbProcLsPoints[i];
                    }
                }

                // Remove duplicate points
                vectorField allPoints(nAllPoints, Zero);
                const scalar tol = 0.001*treeBoundBox(allPointsExt).mag();

                nAllPoints = 0;
                for (const point& pt : allPointsExt)
                {
                    bool duplicate = false;
                    for (label i = 0; i < nAllPoints; ++i)
                    {
                        if (pt.dist(allPoints[i]) < tol)
                        {
                            duplicate = true;
                            break;
                        }
                    }

                    if (!duplicate)
                    {
                        allPoints[nAllPoints++] = pt;
                    }
                }

                allPoints.setSize(nAllPoints);

                if (nAllPoints < 5)
                {
                    FatalErrorInFunction
                        << "There are no enough points for quadrics "
                        << "fitting for a point at processor patch"
                        << abort(FatalError);
                }

                // Transform points
                const vector& origin = points[curPoint];
                const vector axis = normalised(result[curPoint]);
                vector dir(allPoints[0] - points[curPoint]);
                dir.removeCollinear(axis);
                dir.normalise();

                const coordinateSystem cs("cs", origin, axis, dir);

                scalarField W(allPoints.size(), 1.0);

                forAll(allPoints, pI)
                {
                    W[pI] = 1.0/magSqr(allPoints[pI] - points[curPoint]);

                    allPoints[pI] = cs.localPosition(allPoints[pI]);
                }

                scalarRectangularMatrix M
                (
                    allPoints.size(),
                    5,
                    0.0
                );

                for (label i=0; i<allPoints.size(); ++i)
                {
                    M[i][0] = sqr(allPoints[i].x());
                    M[i][1] = sqr(allPoints[i].y());
                    M[i][2] = allPoints[i].x()*allPoints[i].y();
                    M[i][3] = allPoints[i].x();
                    M[i][4] = allPoints[i].y();
                }

                scalarSquareMatrix MtM(5, Zero);

                for (label i = 0; i < MtM.n(); ++i)
                {
                    for (label j = 0; j < MtM.m(); ++j)
                    {
                        for (label k = 0; k < M.n(); ++k)
                        {
                            MtM[i][j] += M[k][i]*M[k][j]*W[k];
                        }
                    }
                }

                scalarField MtR(5, Zero);

                for (label i = 0; i < MtR.size(); ++i)
                {
                    for (label j = 0; j < M.n(); ++j)
                    {
                        MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                    }
                }

                LUsolve(MtM, MtR);

                vector curNormal = vector(MtR[3], MtR[4], -1);

                curNormal = cs.globalVector(curNormal);

                curNormal *= sign(curNormal&result[curPoint]);

                result[curPoint] = curNormal;
            }
        }
    }

    // Correct global points
    if (globalData().nGlobalPoints() > 0)
    {
        const labelList& spLabels = globalData().sharedPointLabels();

        const labelList& addr = globalData().sharedPointAddr();

        for (label k=0; k<globalData().nGlobalPoints(); ++k)
        {
            List<List<vector>> procLsPoints(Pstream::nProcs());

            const label curSharedPointIndex = addr.find(k);

            scalar tol = 0.0;

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                const labelHashSet faceSet(pointFaces[curPoint]);
                const labelList curFaces(faceSet.toc());

                labelHashSet pointSet;
                for (const label facei : curFaces)
                {
                    const labelList& facePoints = faces[facei];
                    pointSet.insert(facePoints);
                }
                pointSet.erase(curPoint);
                labelList curPoints = pointSet.toc();

                vectorField locPoints(points, curPoints);

                procLsPoints[Pstream::myProcNo()] = locPoints;

                tol = 0.001*treeBoundBox(locPoints).mag();
            }

            Pstream::allGatherList(procLsPoints);

            if (curSharedPointIndex != -1)
            {
                label curPoint = spLabels[curSharedPointIndex];

                label nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    nAllPoints += procLsPoints[procI].size();
                }

                vectorField allPoints(nAllPoints, Zero);

                nAllPoints = 0;
                forAll(procLsPoints, procI)
                {
                    forAll(procLsPoints[procI], pointI)
                    {
                        bool duplicate = false;
                        for (label i=0; i<nAllPoints; ++i)
                        {
                            if
                            (
                                mag
                                (
                                    allPoints[i]
                                  - procLsPoints[procI][pointI]
                                )
                              < tol
                            )
                            {
                                duplicate = true;
                                break;
                            }
                        }

                        if (!duplicate)
                        {
                            allPoints[nAllPoints++] =
                                procLsPoints[procI][pointI];
                        }
                    }
                }

                allPoints.setSize(nAllPoints);

                if (nAllPoints < 5)
                {
                    FatalErrorInFunction
                        << "There are no enough points for quadratic "
                        << "fitting for a global processor point "
                        << abort(FatalError);
                }

                // Transform points
                const vector& origin = points[curPoint];
                const vector axis = normalised(result[curPoint]);
                vector dir(allPoints[0] - points[curPoint]);
                dir.removeCollinear(axis);
                dir.normalise();

                coordinateSystem cs("cs", origin, axis, dir);

                scalarField W(allPoints.size(), 1.0);

                forAll(allPoints, pointI)
                {
                    W[pointI]=
                        1.0/magSqr(allPoints[pointI] - points[curPoint]);

                    allPoints[pointI] = cs.localPosition(allPoints[pointI]);
                }

                scalarRectangularMatrix M
                (
                    allPoints.size(),
                    5,
                    0.0
                );

                for (label i=0; i<allPoints.size(); ++i)
                {
                    M[i][0] = sqr(allPoints[i].x());
                    M[i][1] = sqr(allPoints[i].y());
                    M[i][2] = allPoints[i].x()*allPoints[i].y();
                    M[i][3] = allPoints[i].x();
                    M[i][4] = allPoints[i].y();
                }

                scalarSquareMatrix MtM(5, Zero);
                for (label i = 0; i < MtM.n(); ++i)
                {
                    for (label j = 0; j < MtM.m(); ++j)
                    {
                        for (label k = 0; k < M.n(); ++k)
                        {
                            MtM[i][j] += M[k][i]*M[k][j]*W[k];
                        }
                    }
                }

                scalarField MtR(5, Zero);
                for (label i = 0; i < MtR.size(); ++i)
                {
                    for (label j = 0; j < M.n(); ++j)
                    {
                        MtR[i] += M[j][i]*allPoints[j].z()*W[j];
                    }
                }

                LUsolve(MtM, MtR);

                vector curNormal = vector(MtR[3], MtR[4], -1);

                curNormal = cs.globalVector(curNormal);

                curNormal *= sign(curNormal&result[curPoint]);

                result[curPoint] = curNormal;
            }
        }
    }

    for (vector& n : result)
    {
        n.normalise();
    }
}


Foam::tmp<Foam::edgeScalarField> Foam::faMesh::edgeLengthCorrection() const
{
    DebugInFunction
        << "Calculating edge length correction" << endl;

    auto tcorrection = tmp<edgeScalarField>::New
    (
        IOobject
        (
            "edgeLengthCorrection",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless
    );
    auto& correction = tcorrection.ref();

    const vectorField& pointNormals = pointAreaNormals();

    const auto angleCorrection =
        [](const vector& a, const vector& b) -> scalar
        {
            return Foam::cos(0.5*Foam::asin(Foam::mag(a ^ b)));
        };


    // Internal
    {
        auto& fld = correction.primitiveFieldRef();

        forAll(fld, edgei)
        {
            fld[edgei] = angleCorrection
            (
                pointNormals[edges_[edgei].start()],
                pointNormals[edges_[edgei].end()]
            );
        }
    }

    // Boundary
    {
        auto& bfld = correction.boundaryFieldRef();

        forAll(boundary(), patchi)
        {
            const faPatch& fap = boundary()[patchi];
            scalarField& pfld = bfld[patchi];

            label edgei = fap.start();

            forAll(pfld, patchEdgei)
            {
                pfld[patchEdgei] = angleCorrection
                (
                    pointNormals[edges_[edgei].start()],
                    pointNormals[edges_[edgei].end()]
                );

                ++edgei;
            }
        }
    }

    return tcorrection;
}


Foam::tmp<Foam::edgeVectorField> Foam::faMesh::unitLe() const
{
    auto tunitVectors = tmp<edgeVectorField>::New
    (
        IOobject
        (
            "unit(Le)",
            mesh().pointsInstance(),
            faMesh::meshSubDir,
            faMesh::thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        *this,
        dimless,
        (this->Le() / this->magLe())
    );

    // The above is not quite correct when a min-length limiter
    // has been imposed on both Le() and magLe() fields

    // tunitVectors.ref().oriented() = this->Le().oriented();
    return tunitVectors;
}


// ************************************************************************* //
