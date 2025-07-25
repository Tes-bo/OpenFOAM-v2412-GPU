/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "distanceSurface.H"
#include "dictionary.H"
#include "volFields.H"
#include "volPointInterpolation.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(distanceSurface, 0);
}


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::distanceSurface::topologyFilterType
>
Foam::distanceSurface::topoFilterNames_
({
    { topologyFilterType::NONE, "none" },
    { topologyFilterType::LARGEST_REGION, "largestRegion" },
    { topologyFilterType::NEAREST_POINTS, "nearestPoints" },
    { topologyFilterType::PROXIMITY_REGIONS, "proximityRegions" },
    { topologyFilterType::PROXIMITY_FACES, "proximityFaces" },
    { topologyFilterType::PROXIMITY_FACES, "proximity" },
});


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Check that all point hits are valid
static inline void checkAllHits(const UList<pointIndexHit>& nearest)
{
    label notHit = 0;
    for (const pointIndexHit& pHit : nearest)
    {
        if (!pHit.hit())
        {
            ++notHit;
        }
    }

    if (notHit)
    {
        FatalErrorInFunction
            << "Had " << notHit << " faces/cells from "
            << nearest.size() << " without a point hit." << nl
            << "May be caused by a severely degenerate input surface" << nl
            << abort(FatalError);
    }
}


// Normal distance from surface hit point to a point in the mesh
static inline scalar normalDistance_zero
(
    const point& pt,
    const pointIndexHit& pHit,
    const vector& norm
)
{
    const vector diff(pt - pHit.point());

    return (diff & norm);
}


// Signed distance from surface hit point to a point in the mesh,
// the sign is dictated by the normal
static inline scalar normalDistance_nonzero
(
    const point& pt,
    const pointIndexHit& pHit,
    const vector& norm
)
{
    const vector diff(pt - pHit.point());
    const scalar normDist = (diff & norm);

    return Foam::sign(normDist) * Foam::mag(diff);
}


// Normal distance from surface hit point to a point in the mesh
static inline void calcNormalDistance_zero
(
    scalarField& distance,
    const pointField& points,
    const List<pointIndexHit>& nearest,
    const vectorField& normals
)
{
    forAll(nearest, i)
    {
        distance[i] =
            normalDistance_zero(points[i], nearest[i], normals[i]);
    }
}


// Signed distance from surface hit point -> point in the mesh,
// the sign is dictated by the normal
static inline void calcNormalDistance_nonzero
(
    scalarField& distance,
    const pointField& points,
    const List<pointIndexHit>& nearest,
    const vectorField& normals
)
{
    forAll(nearest, i)
    {
        distance[i] =
            normalDistance_nonzero(points[i], nearest[i], normals[i]);
    }
}


// Close to the surface: normal distance from surface hit point
// Far from surface: distance from surface hit point
//
// Note
// This switch may be helpful when working directly with
// distance/gradient fields. Has low overhead otherwise.
// May be replaced in the future (2020-11)
static inline void calcNormalDistance_filtered
(
    scalarField& distance,
    const bitSet& ignoreLocation,
    const pointField& points,
    const List<pointIndexHit>& nearest,
    const vectorField& normals
)
{
    forAll(nearest, i)
    {
        if (ignoreLocation.test(i))
        {
            distance[i] =
                normalDistance_nonzero(points[i], nearest[i], normals[i]);
        }
        else
        {
            distance[i] =
                normalDistance_zero(points[i], nearest[i], normals[i]);
        }
    }
}


// Flat surfaces (eg, a plane) have an extreme change in
// the normal at the edge, which creates a zero-crossing
// extending to infinity.
//
// Ad hoc treatment: require that the surface hit point is within a
// somewhat generous bounding box for the cell (+10%)
//
// Provisioning for filtering based on the cell points,
// but its usefulness remains to be seen (2020-12-09)
template<bool WantPointFilter = false>
static bitSet simpleGeometricFilter
(
    bitSet& ignoreCells,
    const List<pointIndexHit>& nearest,
    const polyMesh& mesh,
    const scalar boundBoxInflate = 0.1  // 10% to catch corners
)
{
    // A deny filter. Initially false (accept everything)
    ignoreCells.resize(mesh.nCells());

    bitSet ignorePoints;
    if (WantPointFilter)
    {
        // Create deny filter
        ignorePoints.resize(mesh.nPoints(), true);
    }

    forAll(nearest, celli)
    {
        const point& pt = nearest[celli].point();

        boundBox cellBb(mesh.cellBb(celli));
        cellBb.inflate(boundBoxInflate);

        if (!cellBb.contains(pt))
        {
            ignoreCells.set(celli);
        }
        else if (WantPointFilter)
        {
            // Good cell candidate, do not ignore its points
            ignorePoints.unset(mesh.cellPoints(celli));
        }
    }

    return ignorePoints;
}


} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::distanceSurface::distanceSurface
(
    const word& defaultSurfaceName,
    const polyMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    geometryPtr_
    (
        searchableSurface::New
        (
            dict.get<word>("surfaceType"),
            IOobject
            (
                dict.getOrDefault("surfaceName", defaultSurfaceName),
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dict
        )
    ),
    distance_(dict.getOrDefault<scalar>("distance", 0)),
    withZeroDistance_(equal(distance_, 0)),
    withSignDistance_
    (
        withZeroDistance_
     || (distance_ < 0)
     || dict.getOrDefault<bool>("signed", true)
    ),

    isoParams_
    (
        dict,
        isoSurfaceParams::ALGO_TOPO,
        isoSurfaceParams::filterType::DIAGCELL
    ),
    topoFilter_
    (
        topoFilterNames_.getOrDefault
        (
            "topology",
            dict,
            topologyFilterType::NONE
        )
    ),
    nearestPoints_(),
    maxDistanceSqr_(Foam::sqr(GREAT)),
    absProximity_(dict.getOrDefault<scalar>("absProximity", 1e-5)),
    cellDistancePtr_(nullptr),
    pointDistance_(),
    surface_(),
    meshCells_(),
    isoSurfacePtr_(nullptr)
{
    if (topologyFilterType::NEAREST_POINTS == topoFilter_)
    {
        dict.readEntry("nearestPoints", nearestPoints_);
    }

    if (dict.readIfPresent("maxDistance", maxDistanceSqr_))
    {
        maxDistanceSqr_ = Foam::sqr(maxDistanceSqr_);
    }
}


Foam::distanceSurface::distanceSurface
(
    const polyMesh& mesh,
    const word& surfaceType,
    const word& surfaceName,
    const isoSurfaceParams& params,
    const bool interpolate
)
:
    distanceSurface
    (
        mesh,
        interpolate,
        surfaceType,
        surfaceName,
        scalar(0),
        true, // redundant - must be signed
        params
    )
{}


Foam::distanceSurface::distanceSurface
(
    const polyMesh& mesh,
    const bool interpolate,
    const word& surfaceType,
    const word& surfaceName,
    const scalar distance,
    const bool useSignedDistance,
    const isoSurfaceParams& params
)
:
    distanceSurface
    (
        mesh,
        interpolate,
        searchableSurface::New
        (
            surfaceType,
            IOobject
            (
                surfaceName,            // name
                mesh.time().constant(), // directory
                "triSurface",           // instance
                mesh.time(),            // registry
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            dictionary()
        ),
        distance,
        useSignedDistance,
        params
    )
{}


Foam::distanceSurface::distanceSurface
(
    const polyMesh& mesh,
    const bool interpolate,
    autoPtr<searchableSurface>&& surface,
    const scalar distance,
    const bool useSignedDistance,
    const isoSurfaceParams& params
)
:
    mesh_(mesh),
    geometryPtr_(surface),
    distance_(distance),
    withZeroDistance_(equal(distance_, 0)),
    withSignDistance_
    (
        withZeroDistance_
     || (distance_ < 0)
     || useSignedDistance
    ),

    isoParams_(params),
    topoFilter_(topologyFilterType::NONE),
    nearestPoints_(),
    maxDistanceSqr_(Foam::sqr(GREAT)),
    absProximity_(1e-5),
    cellDistancePtr_(nullptr),
    pointDistance_(),
    surface_(),
    meshCells_(),
    isoSurfacePtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::distanceSurface::createGeometry()
{
    if (debug)
    {
        Pout<< "distanceSurface::createGeometry updating geometry." << endl;
    }

    // Clear any previously stored topologies
    isoSurfacePtr_.reset(nullptr);
    surface_.clear();
    meshCells_.clear();

    // Doing searches on this surface
    const searchableSurface& geom = geometryPtr_();

    const fvMesh& fvmesh = static_cast<const fvMesh&>(mesh_);

    // Distance to cell centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    cellDistancePtr_.reset
    (
        new volScalarField
        (
            IOobject
            (
                "distanceSurface.cellDistance",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            fvmesh,
            dimensionedScalar(dimLength, GREAT)
        )
    );
    auto& cellDistance = *cellDistancePtr_;


    // For distance = 0 we apply additional geometric filtering
    // to limit the extent of open edges.
    //
    // Does not work with ALGO_POINT

    bitSet ignoreCells, ignoreCellPoints;

    const bool filterCells =
    (
        withZeroDistance_
     && isoParams_.algorithm() != isoSurfaceParams::ALGO_POINT
    );


    // Internal field
    {
        const pointField& cc = fvmesh.C();
        scalarField& fld = cellDistance.primitiveFieldRef();

        List<pointIndexHit> nearest;
        geom.findNearest
        (
            cc,
            // Use initialized field (GREAT) to limit search too
            fld,
            nearest
        );
        checkAllHits(nearest);

        // Geometric pre-filtering when distance == 0

        // NOTE (2021-05-31)
        // Can skip the prefilter if we use proximity-regions filter anyhow
        // but it makes the iso algorithm more expensive and doesn't help
        // unless we start relying on area-based weighting for rejecting regions.

        if (filterCells)
        {
            ignoreCellPoints = simpleGeometricFilter<false>
            (
                ignoreCells,
                nearest,
                fvmesh,

                // Inflate bound box.
                // - To catch corners: approx. 10%
                // - Extra generous for PROXIMITY_REGIONS
                //   (extra weighting for 'bad' faces)
                (
                    topologyFilterType::PROXIMITY_REGIONS == topoFilter_
                  ? 1
                  : 0.1
                )
            );
        }

        if (withSignDistance_)
        {
            vectorField norms;
            geom.getNormal(nearest, norms);

            if (filterCells)
            {
                // With inside/outside switching (see note above)
                calcNormalDistance_filtered
                (
                    fld,
                    ignoreCells,
                    cc,
                    nearest,
                    norms
                );
            }
            else if (withZeroDistance_)
            {
                calcNormalDistance_zero(fld, cc, nearest, norms);
            }
            else
            {
                calcNormalDistance_nonzero(fld, cc, nearest, norms);
            }
        }
        else
        {
            calcAbsoluteDistance(fld, cc, nearest);
        }
    }


    // Patch fields
    {
        forAll(fvmesh.C().boundaryField(), patchi)
        {
            const pointField& cc = fvmesh.C().boundaryField()[patchi];
            scalarField& fld = cellDistance.boundaryFieldRef()[patchi];

            List<pointIndexHit> nearest;
            geom.findNearest
            (
                cc,
                scalarField(cc.size(), GREAT),
                nearest
            );
            checkAllHits(nearest);

            if (withSignDistance_)
            {
                vectorField norms;
                geom.getNormal(nearest, norms);

                if (withZeroDistance_)
                {
                    // Slight inconsistency in boundary vs interior when
                    // cells are filtered, but the patch fields are only
                    // used by isoSurfacePoint, and filtering is disabled
                    // for that anyhow.

                    calcNormalDistance_zero(fld, cc, nearest, norms);
                }
                else
                {
                    calcNormalDistance_nonzero(fld, cc, nearest, norms);
                }
            }
            else
            {
                calcAbsoluteDistance(fld, cc, nearest);
            }
        }
    }


    // On processor patches the mesh.C() will already be the cell centre
    // on the opposite side so no need to swap cellDistance.


    // Distance to points
    pointDistance_.resize(fvmesh.nPoints());
    pointDistance_ = GREAT;
    {
        const pointField& pts = fvmesh.points();
        scalarField& fld = pointDistance_;

        List<pointIndexHit> nearest;
        geom.findNearest
        (
            pts,
            // Use initialized field (GREAT) to limit search too
            pointDistance_,
            nearest
        );
        checkAllHits(nearest);

        if (withSignDistance_)
        {
            vectorField norms;
            geom.getNormal(nearest, norms);

            if (filterCells)
            {
                // With inside/outside switching (see note above)
                calcNormalDistance_filtered
                (
                    fld,
                    ignoreCellPoints,
                    pts,
                    nearest,
                    norms
                );
            }
            else if (withZeroDistance_)
            {
                calcNormalDistance_zero(fld, pts, nearest, norms);
            }
            else
            {
                calcNormalDistance_nonzero(fld, pts, nearest, norms);
            }
        }
        else
        {
            calcAbsoluteDistance(fld, pts, nearest);
        }
    }


    // Don't need ignoreCells if there is nothing to ignore.
    if (ignoreCells.none())
    {
        ignoreCells.clearStorage();
    }
    else if (filterCells && topologyFilterType::NONE != topoFilter_)
    {
        // For refine blocked cells (eg, checking actual cells cut)
        isoSurfaceBase isoCutter
        (
            mesh_,
            cellDistance,
            pointDistance_,
            distance_
        );

        if (topologyFilterType::LARGEST_REGION == topoFilter_)
        {
            refineBlockedCells(ignoreCells, isoCutter);
            filterKeepLargestRegion(ignoreCells);
        }
        else if (topologyFilterType::NEAREST_POINTS == topoFilter_)
        {
            refineBlockedCells(ignoreCells, isoCutter);
            filterKeepNearestRegions(ignoreCells);
        }

        // Note: apply similar filtering for PROXIMITY_REGIONS later instead
    }

    // Don't need point filter beyond this point
    ignoreCellPoints.clearStorage();


    if (debug)
    {
        Pout<< "Writing cell distance:" << cellDistance.objectPath() << endl;
        cellDistance.write();

        pointScalarField pDist
        (
            IOobject
            (
                "distanceSurface.pointDistance",
                fvmesh.time().timeName(),
                fvmesh.time(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            pointMesh::New(fvmesh),
            dimensionedScalar(dimLength, Zero)
        );
        pDist.primitiveFieldRef() = pointDistance_;

        Pout<< "Writing point distance:" << pDist.objectPath() << endl;
        pDist.write();
    }

    isoSurfacePtr_.reset
    (
        isoSurfaceBase::New
        (
            isoParams_,
            cellDistance,
            pointDistance_,
            distance_,
            ignoreCells
        )
    );


    // Restrict ignored cells to those actually cut
    if (filterCells && topologyFilterType::PROXIMITY_REGIONS == topoFilter_)
    {
        isoSurfaceBase isoCutter
        (
            mesh_,
            cellDistance,
            pointDistance_,
            distance_
        );

        refineBlockedCells(ignoreCells, isoCutter);
    }


    // ALGO_POINT still needs cell, point fields (for interpolate)
    // The others can do straight transfer

    // But also flatten into a straight transfer for proximity filtering
    if
    (
        isoParams_.algorithm() != isoSurfaceParams::ALGO_POINT
     || topologyFilterType::PROXIMITY_FACES == topoFilter_
     || topologyFilterType::PROXIMITY_REGIONS == topoFilter_
    )
    {
        surface_.transfer(static_cast<meshedSurface&>(*isoSurfacePtr_));
        meshCells_.transfer(isoSurfacePtr_->meshCells());

        isoSurfacePtr_.reset(nullptr);
        cellDistancePtr_.reset(nullptr);
        pointDistance_.clear();
    }

    if (topologyFilterType::PROXIMITY_FACES == topoFilter_)
    {
        filterFaceProximity();
    }
    else if (topologyFilterType::PROXIMITY_REGIONS == topoFilter_)
    {
        filterRegionProximity(ignoreCells);
    }

    if (debug)
    {
        print(Pout, debug);
        Pout<< endl;
    }
}


void Foam::distanceSurface::print(Ostream& os, int level) const
{
    os  << " surface:" << surfaceName()
        << " distance:" << distance()
        << " topology:" << topoFilterNames_[topoFilter_];

    isoParams_.print(os);

    if (level)
    {
        os  << "  faces:" << surface().surfFaces().size()
            << "  points:" << surface().points().size();
    }
}


// ************************************************************************* //
