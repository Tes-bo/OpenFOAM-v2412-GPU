/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::searchableSurfaces

Description
    Container for searchableSurfaces.
    The collection is specified as a dictionary. For example,

    \verbatim
    geometry
    {
        surface1
        {
            type  ...;
        }
        surface2
        {
            type  ...;
        }
    }
    \endverbatim

    \heading The Sub-dictionary parameters
    \table
        Property    | Description                       | Required  | Default
        name        | alternative name for surface      | no    | dict name
        regions     | Region names sub-dictionary       | no    |
    \endtable

SourceFiles
    searchableSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_searchableSurfaces_H
#define Foam_searchableSurfaces_H

#include "searchableSurface.H"
#include "labelPair.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class coordSetWriter;
class triSurface;

/*---------------------------------------------------------------------------*\
                           Class searchableSurfaces Declaration
\*---------------------------------------------------------------------------*/

class searchableSurfaces
:
    public PtrList<searchableSurface>
{
    // Private data

        //- Surface names
        wordList names_;

        //- Region names per surface
        List<wordList> regionNames_;

        //- Indices of all surfaces. Precalculated and stored.
        labelList allSurfaces_;


    // Private Member Functions

        //- Is edge on face
        static bool connected
        (
            const triSurface& s,
            const label edgeI,
            const pointIndexHit& hit
        );


        //- No copy construct
        searchableSurfaces(const searchableSurfaces&) = delete;

        //- No copy assignment
        void operator=(const searchableSurfaces&) = delete;


public:

    ClassName("searchableSurfaces");


    // Constructors

        //- Construct with length specified. Fill later.
        explicit searchableSurfaces(const label);


        ////- Construct from list of dictionaries
        //searchableSurfaces(const IOobject&, const PtrList<dictionary>&);

        //- Construct from dictionary.
        //  \param singleRegionName controls if names are constructed as
        //      surfaceName "_" regionName (singleRegionName false) or for
        //      single region surfaces as surfaceName only (singleRegionName
        //      true)
        searchableSurfaces
        (
            const IOobject&,
            const dictionary&,
            const bool singleRegionName
        );


    // Member Functions

        //- Surface names, not region names
        const wordList& names() const
        {
            return names_;
        }

        //- Surface names, not region names
        wordList& names()
        {
            return names_;
        }

        //- Region names per surface
        const List<wordList>& regionNames() const
        {
            return regionNames_;
        }

        //- Region names per surface
        List<wordList>& regionNames()
        {
            return regionNames_;
        }


        ////- If adding surfaces 'by hand'
        //HashTable<labelPair>& regionNames()
        //{
        //    return regionNames_;
        //}
        ////- Get surface and region for a name
        //const labelPair& surfaceRegion(const word& globalRegion) const
        //{
        //    return regionNames_[globalRegion];
        //}

        //- Find index of surface. Return -1 if not found.
        label findSurfaceID(const word& name) const;

        label findSurfaceRegionID
        (
            const word& surfaceName,
            const word& regionName
        ) const;

        // Multiple point queries.

            //- Find any intersection. Return hit point information and
            //  surface number. If multiple surfaces hit the first surface
            //  is returned, not necessarily the nearest (to start).
            void findAnyIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            //- Find all intersections in order from start to end. Returns for
            //  every hit the surface and the hit info.
            void findAllIntersections
            (
                const pointField& start,
                const pointField& end,
                labelListList& surfaces,
                List<List<pointIndexHit>>&
            ) const;

            //Find intersections of edge nearest to both endpoints.
            void findNearestIntersection
            (
                const pointField& start,
                const pointField& end,
                labelList& surface1,
                List<pointIndexHit>& hit1,
                labelList& surface2,
                List<pointIndexHit>& hit2
            ) const;

            //- Find nearest. Return -1 (and a miss()) or surface and nearest
            //  point.
            void findNearest
            (
                const pointField&,
                const scalarField& nearestDistSqr,
                labelList& surfaces,
                List<pointIndexHit>&
            ) const;

            void findNearest
            (
                const labelListList& regionIndices,
                const pointField& samples,
                const scalarField& nearestDistSqr,
                labelList& nearestSurfaces,
                List<pointIndexHit>& nearestInfo
            ) const;

            //- Calculate bounding box
            boundBox bounds() const;

        // Checking

            //- Are all surfaces closed and manifold
            bool checkClosed(const bool report) const;

            //- Are all (triangulated) surfaces consistent normal orientation
            bool checkNormalOrientation(const bool report) const;

            //- Are all bounding boxes of similar size
            bool checkSizes(const scalar maxRatio, const bool report) const;

            //- Do surfaces self-intersect or intersect others
            bool checkIntersection
            (
                const scalar tol,
                autoPtr<coordSetWriter>& setWriter,
                const bool report
            ) const;

            //- Check triangle quality
            bool checkQuality
            (
                const scalar minQuality,
                const bool report
            ) const;

            //- All topological checks. Return number of failed checks
            label checkTopology(const bool report) const;

            //- All geometric checks. Return number of failed checks
            label checkGeometry
            (
                const scalar maxRatio,
                const scalar tolerance,
                autoPtr<coordSetWriter>& setWriter,
                const scalar minQuality,
                const bool report
            ) const;

            //- Write some stats
            void writeStats(const List<wordList>&, Ostream&) const;


    // Member Operators

        //- Return const and non-const reference to searchableSurface by index.
        using PtrList<searchableSurface>::operator[];

        //- Return const reference to searchableSurface by name.
        const searchableSurface& operator[](const word&) const;

        //- Return reference to searchableSurface by name.
        searchableSurface& operator[](const word&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
