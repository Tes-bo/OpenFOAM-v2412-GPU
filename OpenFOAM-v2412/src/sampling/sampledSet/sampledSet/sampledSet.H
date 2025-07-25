/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::sampledSet

Group
    grpUtilitiesFunctionObjects

Description
    Holds list of sampling points which is filled at construction time.
    Various implementations of this base class to e.g. get sampling points
    at uniform distance along a line (uniformSet) or directly specified
    (cloudSet)

    Each 'sampledSet' has a name and a specifier of how the axis should be
    write (x/y/z component or all 3 components)

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        axis     | x, y, z, xyz, distance                  | yes      |
    \endtable

SourceFiles
    sampledSet.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_sampledSet_H
#define Foam_sampledSet_H

#include "coordSet.H"
#include "autoPtr.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class meshSearch;

/*---------------------------------------------------------------------------*\
                         Class sampledSet Declaration
\*---------------------------------------------------------------------------*/

class sampledSet
:
    public coordSet
{
    // Private Data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Reference to mesh searching class
        const meshSearch& searchEngine_;


protected:

    // Protected Data

        //- Segment numbers
        labelList segments_;

        //- Cell numbers
        labelList cells_;

        //- Face numbers (-1 if not known)
        labelList faces_;


    // Protected Member Functions

        //- Check for consistent sizing
        void checkDimensions() const;

        //- Returns cell next to boundary face
        label getBoundaryCell(const label) const;

        //- Returns the neighbour cell or the owner if face in on the boundary
        label getNeighbourCell(const label) const;

        //- Return the cell in which the point on the sample line
        //  resides if found otherwise return -1
        label pointInCell(const point& p, const label samplei) const;

        //- Calculates inproduct of face normal and vector sample-face centre
        //  <0 if sample inside.
        scalar calcSign(const label facei, const point& sample) const;

        //- Returns face label (or -1) of face which is close to sample
        label findNearFace
        (
            const label celli,
            const point& sample,
            const scalar smallDist
        ) const;

        //- Moves sample in direction of -n to it is 'inside' of facei
        point pushIn
        (
            const point& sample,
            const label facei
        ) const;

        //- Calculates start of tracking given samplePt and first boundary
        //  intersection
        //  (bPoint, bFacei) (bFacei == -1 if no boundary intersection)
        //  Returns true if trackPt is valid sampling point. Sets trackPt,
        //  trackFacei, trackCelli (-1 if no tracking point found)
        bool getTrackingPoint
        (
            const point& samplePt,
            const point& bPoint,
            const label bFacei,
            const scalar smallDist,

            point& trackPt,
            label& trackCelli,
            label& trackFacei
        ) const;

        //- Set sample data. Copy list contents.
        void setSamples
        (
            const List<point>& samplingPts,
            const labelList& samplingCells,
            const labelList& samplingFaces,
            const labelList& samplingSegments,
            const scalarList& samplingDistance
        );

        //- Set sample data. Move list contents.
        void setSamples
        (
            List<point>&& samplingPts,
            labelList&& samplingCells,
            labelList&& samplingFaces,
            labelList&& samplingSegments,
            scalarList&& samplingDistance
        );

public:

    //- Runtime type information
    TypeName("sampledSet");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sampledSet,
            word,
            (
                const word& name,
                const polyMesh& mesh,
                const meshSearch& searchEngine,
                const dictionary& dict
            ),
            (name, mesh, searchEngine, dict)
        );


        //- Class used for the read-construction of
        //  PtrLists of sampledSet
        class iNew
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;
            const meshSearch& search_;

        public:

            iNew(const polyMesh& mesh, const meshSearch& searchEngine)
            :
                mesh_(mesh),
                search_(searchEngine)
            {}

            autoPtr<sampledSet> operator()(Istream& is) const
            {
                word name(is);
                dictionary dict(is);
                return sampledSet::New(name, mesh_, search_, dict);
            }
        };


        //- PtrList read-construction helper that captures dictionaries used
        //- during creation.
        class iNewCapture
        {
            //- Reference to the volume mesh
            const polyMesh& mesh_;
            const meshSearch& search_;

            //- Captured (recorded) dictionaries
            DynamicList<dictionary>& capture_;

        public:

            iNewCapture
            (
                const polyMesh& mesh,
                const meshSearch& searchEngine,
                DynamicList<dictionary>& capture
            )
            :
                mesh_(mesh),
                search_(searchEngine),
                capture_(capture)
            {}

            autoPtr<sampledSet> operator()(Istream& is) const
            {
                word name(is);
                dictionary& dict = capture_.emplace_back(is);

                return sampledSet::New(name, mesh_, search_, dict);
            }
        };


    // Constructors

        //- Construct from components
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const coordSet::coordFormat axisType
        );

        //- Construct from components
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis
        );

        //- Construct from dictionary
        sampledSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );

        //- Clone
        autoPtr<sampledSet> clone() const
        {
            NotImplemented;
            return nullptr;
        }


    // Selectors

        //- Return a reference to the selected sampledSet
        static autoPtr<sampledSet> New
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~sampledSet() = default;


    // Member Functions

        const polyMesh& mesh() const noexcept
        {
            return mesh_;
        }

        const meshSearch& searchEngine() const noexcept
        {
            return searchEngine_;
        }

        const labelList& segments() const noexcept
        {
            return segments_;
        }

        const labelList& cells() const noexcept
        {
            return cells_;
        }

        const labelList& faces() const noexcept
        {
            return faces_;
        }

        //- Output for debugging
        Ostream& write(Ostream&) const;


    // Other

        /// //- Gather and sort.
        /// \return (on master) gathered set and overall sort order
        /// autoPtr<coordSet> gather
        /// (
        ///    labelList& sortOrder,
        ///    labelList& allSegments
        /// ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
