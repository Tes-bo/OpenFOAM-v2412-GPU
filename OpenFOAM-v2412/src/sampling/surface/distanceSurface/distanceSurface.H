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

Class
    Foam::distanceSurface

Description
    A surface defined by a distance from an input searchable surface.
    Uses an iso-surface algorithm (cell, topo, point) for constructing the
    distance surface.

    For a zero-distance surface, it performs additional checks and supports
    filtering to handle the surface boundaries.

Usage
    Example of function object partial specification:
    \verbatim
    surfaces
    {
        surface1
        {
            type        distanceSurface;
            surfaceType triSurfaceMesh;
            surfaceName something.obj;
            topology    proximityFaces;
        }

        surface2
        {
            type        distanceSurface;
            surfaceType triSurfaceMesh;
            surfaceName other.obj;

            topology    nearestPoints;
            nearestPoints
            (
                (0 0 0)
                (10 10 0)
            );

            // Max search distance for nearestPoints
            maxDistance 0.005;
        }
    }
    \endverbatim

    Dictionary controls:
    \table
        Property | Description                              | Required | Default
        distance | distance from surface                    | no  | 0
        signed   | Use sign when distance is positive       | no  | true
        isoMethod | Iso-algorithm (cell/topo/point)         | no  | default
        regularise | Face simplification (enum or bool)     | no  | true
        bounds   | Limit with bounding box                  | no  |
        surfaceType | Type of surface                       | yes |
        surfaceName | Name of surface in \c triSurface/     | no  | dict name
        topology    | Topology filter name                  | no  | none
        nearestPoints | Points for point-based segmentation | no  |
        maxDistance | Max search distance for nearestPoints | no  | GREAT
        absProximity | Max proximity of face centres        | no  | 1e-5
    \endtable

    Topology/Filtering (for zero-distance only).
    These represent different ways to tackle the "ragged edge" problem.

    - \c none :
        No filtering

    - \c proximityFaces or \c proximity (post-filter):
        Checks the resulting faces against the original search surface
        and rejects faces with a distance greater than \c absProximity.

    - \c proximityRegions (post-filter):
        Checks the distance of the resulting faces against the original
        search surface. Filters based on the area-weighted distance
        of each topologically connected region.
        If the area-weighted distance of a region is greater than
        \c absProximity, the entire region is rejected.

    - \c largestRegion (pre-filter):
        The cut cells are checked for topological connectivity and the
        region with the most number of cut cells is retained.

    - \c nearestPoints (pre-filter):
        The cut cells split into regions, the regions closest to the
        user-defined points are retained.
        Uses \c maxDistance for additional control.
    .

Note
    For distance = 0, some special adjustments.
    - Always signed (ignoring the input value).
    - Use normal distance from surface (for better treatment of open edges).
    - Additional checks for open surfaces edges are used to limit the extend
      of resulting distance surface.
      The resulting surface elements will, however, contain partial cell
      coverage. NB: Not applicable if the \c point isoMethod is used.

The keyword \c cell (bool value) which was use in 1906 and earlier to switch
between point/cell algorithms is now ignored (2020-12).

Changed default algorithm from cell to topo (2020-12).

SourceFiles
    distanceSurface.C
    distanceSurfaceFilter.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_distanceSurface_H
#define Foam_distanceSurface_H

#include "sampledSurface.H"
#include "searchableSurface.H"
#include "isoSurfaceBase.H"
#include "pointIndexHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class distanceSurface Declaration
\*---------------------------------------------------------------------------*/

class distanceSurface
{
    // Data Types

        //- The type of pre/post face-filtering
        enum class topologyFilterType : uint8_t
        {
            NONE,               //!< No additional filtering
            LARGEST_REGION,     //!< Retain largest region
            NEAREST_POINTS,     //!< Retain regions nearest to the points
            PROXIMITY_REGIONS,  //!< Retain regions with good surface proximity
            PROXIMITY_FACES,    //!< Retain faces with good surface proximity
            PROXIMITY = PROXIMITY_FACES
        };

        //- Names for the topology filter
        static const Enum<topologyFilterType> topoFilterNames_;


    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Searchable surface
        const autoPtr<searchableSurface> geometryPtr_;

        //- Distance value
        const scalar distance_;

        //- Distance is zero. Implies signed and additional optimizations
        const bool withZeroDistance_;

        //- Use signed distance
        const bool withSignDistance_;

        //- Parameters for iso-surface (algorithm, filter, mergeTol, etc)
        isoSurfaceParams isoParams_;

        //- Optional topology face-filtering
        topologyFilterType topoFilter_;

        //- Points for nearest-points segmentation
        pointField nearestPoints_;

        //- Max search distance squared (for nearestPoints)
        scalar maxDistanceSqr_;

        //- Max distance for proximity check (post-filtering)
        scalar absProximity_;

        //- Distance to cell centres
        autoPtr<volScalarField> cellDistancePtr_;

        //- Distance to points
        scalarField pointDistance_;


    // Sampling geometry. Directly stored or via an iso-surface (ALGO_POINT)

        //- The extracted surface (direct storage)
        mutable meshedSurface surface_;

        //- For every face the original cell in mesh (direct storage)
        mutable labelList meshCells_;

        //- Extracted iso-surface, for interpolators
        mutable autoPtr<isoSurfaceBase> isoSurfacePtr_;


    // Private Member Functions

        //- Absolute distances from hit points
        //  Hit/miss checks have been done elsewhere.
        static inline void calcAbsoluteDistance
        (
            scalarField& distance,
            const pointField& points,
            const List<pointIndexHit>& nearest
        )
        {
            forAll(nearest, i)
            {
                distance[i] = Foam::mag(points[i] - nearest[i].point());
            }
        }


protected:

    // Protected Member Functions

        //- Is currently backed by an isoSurfacePtr_
        bool hasIsoSurface() const
        {
            return bool(isoSurfacePtr_);
        }

        //- Interpolate volume field onto surface points
        template<class Type>
        tmp<Field<Type>> isoSurfaceInterpolate
        (
            const VolumeField<Type>& cellValues,
            const Field<Type>& pointValues
        ) const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->interpolate(cellValues, pointValues);
            }

            return nullptr;
        }


    // Private Member Functions

        //- Re-filter the blocked cells based on the anticipated cuts
        //  Uses a lightweight variant of cutting.
        bool refineBlockedCells
        (
            bitSet& ignoreCells,
            const isoSurfaceBase& isoContext
        ) const;

        //- Prepare blockedFaces for region split
        bitSet filterPrepareRegionSplit(const bitSet& ignoreCells) const;

        //- Keep region with the most cuts (after region split)
        void filterKeepLargestRegion(bitSet& ignoreCells) const;

        //- Keep region(s) closest to the nearest points
        void filterKeepNearestRegions(bitSet& ignoreCells) const;

        //- Remove region(s) that have far faces
        void filterRegionProximity(bitSet& ignoreCells) const;

        //- Adjust extracted iso-surface to remove far faces
        void filterFaceProximity();


public:

    //- Runtime type information
    TypeName("distanceSurface");


    // Constructors

        //- Construct from dictionary
        distanceSurface
        (
            const word& defaultSurfaceName,
            const polyMesh& mesh,
            const dictionary& dict
        );

        //- Construct from components with zero-distanced
        distanceSurface
        (
            const polyMesh& mesh,
            const word& surfaceType,
            const word& surfaceName,
            const isoSurfaceParams& params = isoSurfaceParams(),
            const bool interpolate = false
        );

        //- Construct from components
        distanceSurface
        (
            const polyMesh& mesh,
            const bool interpolate,
            const word& surfaceType,
            const word& surfaceName,
            const scalar distance,
            const bool useSignedDistance,
            const isoSurfaceParams& params = isoSurfaceParams()
        );

        distanceSurface
        (
            const polyMesh& mesh,
            const bool interpolate,
            autoPtr<searchableSurface>&& surface,
            const scalar distance,
            const bool useSignedDistance,
            const isoSurfaceParams& params = isoSurfaceParams()
        );


    //- Destructor
    virtual ~distanceSurface() = default;


    // Member Functions

        //- Create/recreate the distance surface
        void createGeometry();

        //- The name of the underlying searchableSurface
        const word& surfaceName() const
        {
            return geometryPtr_->name();
        }

        //- The distance to the underlying searchableSurface
        scalar distance() const noexcept
        {
            return distance_;
        }

        //- The underlying surface
        const meshedSurface& surface() const
        {
            if (isoSurfacePtr_)
            {
                return *isoSurfacePtr_;
            }
            return surface_;
        }

        //- The underlying surface
        meshedSurface& surface()
        {
            if (isoSurfacePtr_)
            {
                return *isoSurfacePtr_;
            }
            return surface_;
        }

        //- For each face, the original cell in mesh
        const labelList& meshCells() const
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return meshCells_;
        }

        //- For each face, the original cell in mesh
        labelList& meshCells()
        {
            if (isoSurfacePtr_)
            {
                return isoSurfacePtr_->meshCells();
            }
            return meshCells_;
        }


    // Output

        //- Print information
        void print(Ostream& os, int level=0) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
