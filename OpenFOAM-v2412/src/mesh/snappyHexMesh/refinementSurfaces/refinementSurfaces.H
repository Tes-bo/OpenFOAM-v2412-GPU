/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::refinementSurfaces

Description
    Container for data on surfaces used for surface-driven refinement.
    Contains all the data about the level of refinement needed per
    surface.

SourceFiles
    refinementSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef refinementSurfaces_H
#define refinementSurfaces_H

#include "triSurfaceGeoMesh.H"
#include "triSurfaceFields.H"
#include "vectorList.H"
#include "pointIndexHit.H"
#include "surfaceZonesInfo.H"
#include "volumeType.H"
#include "pointList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;
class shellSurfaces;
class triSurfaceMesh;

/*---------------------------------------------------------------------------*\
                           Class refinementSurfaces Declaration
\*---------------------------------------------------------------------------*/

class refinementSurfaces
{
    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are refinement ones
        labelList surfaces_;

        //- Surface name (word)
        wordList names_;

        //- List of surface zone (face and cell zone) information
        PtrList<surfaceZonesInfo> surfZones_;

        //- From surface to starting global region
        labelList regionOffset_;

        //- From global region number to surface
        labelList regionToSurface_;

        //- From global region number to refinement level
        labelList minLevel_;

        //- From global region number to refinement level
        labelList maxLevel_;

        //- From global region number to small-gap level
        labelList gapLevel_;

        //- From global region number to cell level at which blockage detection
        //  needs to apply
        labelList blockLevel_;

        //- From global region number to cell level at which leakage detection
        //  needs to apply
        labelList leakLevel_;

        //- From global region number to small-gap level specification
        List<FixedList<label, 3>> extendedGapLevel_;

        //- From global region number to side of surface to detect
        List<volumeType> extendedGapMode_;

        //- From global region number to whether to allow selfProximity
        //  (in gap refinement)
        boolList selfProximity_;

        //- From global region number to curvature specification
        List<FixedList<label, 4>> extendedCurvatureLevel_;

        //- From global region number to perpendicular angle
        scalarField perpendicularAngle_;

        //- From global region number to patchType
        PtrList<dictionary> patchInfo_;

        //- From global region number to whether to add buffer layers
        boolList addBufferLayers_;

        //- Are we operating in test mode?
        const bool dryRun_;


    // Private Member Functions

        //- Given intersection results with geom detect local shell refinement
        //  level (possibly cached on triangles of geom)
        labelList findHigherLevel
        (
            const searchableSurface& geom,
            const shellSurfaces& shells,
            const List<pointIndexHit>& intersectionInfo,
            const labelList& surfaceLevel
        ) const;

        //- Calculate global region to surface
        static labelList calcSurfaceIndex
        (
            const searchableSurfaces& allGeometry,
            const labelList& surfaces
        );

        //- No copy construct
        refinementSurfaces(const refinementSurfaces&) = delete;

        //- No copy assignment
        void operator=(const refinementSurfaces&) = delete;


public:

    // Constructors

        //- Construct from surfaces and dictionary
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary&,
            const label gapLevelIncrement,
            const bool dryRun
        );

        //- Construct from components
        refinementSurfaces
        (
            const searchableSurfaces& allGeometry,
            const labelList& surfaces,
            const wordList& names,
            const PtrList<surfaceZonesInfo>& surfZones,
            const labelList& regionOffset,
            const labelList& minLevel,
            const labelList& maxLevel,
            const labelList& gapLevel,
            const scalarField& perpendicularAngle,
            PtrList<dictionary>& patchInfo,
            const bool dryRun
        );


    // Member Functions

        // Access

            const searchableSurfaces& geometry() const
            {
                return allGeometry_;
            }

            const labelList& surfaces() const
            {
                return surfaces_;
            }

            //- Names of surfaces
            const wordList& names() const
            {
                return names_;
            }

            const PtrList<surfaceZonesInfo>& surfZones() const
            {
                return surfZones_;
            }

            //- From surface to starting global region
            const labelList& regionOffset() const
            {
                return regionOffset_;
            }

            //- From global region number to refinement level
            const labelList& minLevel() const
            {
                return minLevel_;
            }

            //- From global region number to refinement level
            const labelList& maxLevel() const
            {
                return maxLevel_;
            }

            //- From global region number to small gap refinement level
            const labelList& gapLevel() const
            {
                return gapLevel_;
            }

            //- From global region number to cell level at which blockage
            //- detection is applied
            const labelList& blockLevel() const
            {
                return blockLevel_;
            }

            //- From global region number to cell level at which leakage
            //- detection is applied. labelMax if not set.
            const labelList& leakLevel() const
            {
                return leakLevel_;
            }

            //- From global region number to specification of gap and its
            //  refinement: 3 labels specifying
            //  - minimum wanted number of cells in the gap
            //  - minimum cell level when to start trying to detect gaps
            //  - maximum cell level to refine to (so do not detect gaps if
            //    cell >= maximum level)
            const List<FixedList<label, 3>>& extendedGapLevel() const
            {
                return extendedGapLevel_;
            }

            //- From global region number to side of surface to detect
            const List<volumeType>& extendedGapMode() const
            {
                return extendedGapMode_;
            }

            //- From global region number to whether to detect gaps to same
            //  surface (in gap refinement)
            const boolList& gapSelf() const
            {
                return selfProximity_;
            }

            //- From global region number to specification of curvature
            //  refinement: 4 labels specifying
            //  - minimum wanted number of cells in the curvature radius
            //  - ?minimum cell level when to start trying to detect gaps
            //  - maximum cell level to refine to (so do not detect curvature if
            //    celllevel >= maximum level)
            //  - minimum radius to ignore (expressed as refinement level).
            //    This can be used to ignore feature-edges. Set to -1 to ignore.
            const List<FixedList<label, 4>>& extendedCurvatureLevel() const
            {
                return extendedCurvatureLevel_;
            }

            //- From global region number to perpendicular angle
            const scalarField& perpendicularAngle() const
            {
                return perpendicularAngle_;
            }

            //- From global region number to patch type
            const PtrList<dictionary>& patchInfo() const
            {
                return patchInfo_;
            }

            //- From global region number to whether to add buffer layers
            //  when snapping
            const boolList& addBufferLayers() const
            {
                return addBufferLayers_;
            }


        // Helper

            //- From surface and region on surface to global region
            label globalRegion(const label surfI, const label regionI) const
            {
                return regionOffset_[surfI]+regionI;
            }

            //- From global region to surface + region
            labelPair whichSurface(const label globalRegionI) const;

            //- Min level for surface and region on surface
            label minLevel(const label surfI, const label regionI) const
            {
                return minLevel_[globalRegion(surfI, regionI)];
            }

            //- Max level for surface and region on surface
            label maxLevel(const label surfI, const label regionI) const
            {
                return maxLevel_[globalRegion(surfI, regionI)];
            }

            label nRegions() const
            {
                return minLevel_.size();
            }

            //- Per surface the maximum extendedGapLevel over all its regions
            labelList maxGapLevel() const;

            //- Per surface the maximum curvatureLevel over all its regions
            labelList maxCurvatureLevel() const;

            //- Calculate minLevelFields according to both surface- and
            //  shell-based levels
            void setMinLevelFields
            (
                const shellSurfaces& shells
            );

            //- Update minLevelFields according to (triSurface-only) curvature
            void setCurvatureMinLevelFields
            (
                const scalar cosAngle,
                const scalar level0EdgeLength
            );

            ////- Helper: count number of triangles per region
            //static labelList countRegions(const triSurface&);


        // Searching

            //- Find intersection of edge. Return -1 or first surface
            //  with higher (than currentLevel) minlevel.
            //  Return surface number and level.
            void findHigherIntersection
            (
                const shellSurfaces& shells,

                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                labelList& surfaces,
                labelList& surfaceLevel
            ) const;

            //- Find all intersections of edge with any surface with applicable
            //  min/max refinement level. Unsorted order.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                const labelList& globalMinLevel,    // level per surfregion
                const labelList& globalMaxLevel,    // level per surfregion

                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find all intersections of edge with any surface with applicable
            //  min/max refinement level. Unsorted order.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                const labelList& currentLevel,  // current cell refinement level

                const labelList& globalMinLevel,    // level per surfregion
                const labelList& globalMaxLevel,    // level per surfregion

                List<pointList>& surfaceLocation,
                List<vectorList>& surfaceNormal,
                labelListList& surfaceLevel
            ) const;

            //- Find intersection nearest to the endpoints. surface1,2 are
            //  not indices into surfacesToTest but refinement surface indices.
            //  Returns surface, region on surface (so not global surface)
            //  and position on surface.
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2
            ) const;

            //- findNearestIntersection but also get normals
            void findNearestIntersection
            (
                const labelList& surfacesToTest,
                const pointField& start,
                const pointField& end,

                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& region1,
                vectorField& normal1,

                labelList& surface2,
                List<pointIndexHit>& hit2,
                labelList& region2,
                vectorField& normal2
            ) const;

            //- Find nearest (to start only) intersection of edge
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                vectorField& normal
            ) const;

            //- Find nearest (to start only) intersection of edge
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                labelList& regions,
                List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Find nearest (to start only) intersection of edge
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&,
                vectorField& normal
            ) const;

            //- Used for debugging only: find intersection of edge.
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces.
            void findNearest
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find nearest point on surfaces. Return surface and region on
            //  surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                labelList& hitRegion
            ) const;

            //- Find nearest point on surfaces. Return surface, region and
            //  normal on surface (so not global surface)
            void findNearestRegion
            (
                const labelList& surfacesToTest,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& hitSurface,
                List<pointIndexHit>& hitInfo,
                labelList& hitRegion,
                vectorField& hitNormal
            ) const;

            //- Detect if a point is 'inside' (closed) surfaces.
            //  Returns -1 if not, returns first surface it is.
            void findInside
            (
                const labelList& surfacesToTest,
                const pointField& pt,
                labelList& insideSurfaces
            ) const;

            // Region wise searching

                //- Find nearest point on selected regions of surfaces.
                void findNearest
                (
                    const labelList& surfacesToTest,
                    const labelListList& regions,

                    const pointField& samples,
                    const scalarField& nearestDistSqr,

                    labelList& hitSurface,
                    List<pointIndexHit>& hitInfo
                ) const;

                //- Find nearest point on selected regions of surfaces.
                void findNearestRegion
                (
                    const labelList& surfacesToTest,
                    const labelListList& regions,

                    const pointField& samples,
                    const scalarField& nearestDistSqr,

                    labelList& hitSurface,
                    List<pointIndexHit>& hitInfo,
                    labelList& hitRegion,
                    vectorField& hitNormal
                ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
