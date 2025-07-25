/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::shortestPathSet

Description
    Finds shortest path (in terms of cell centres) to walk on mesh from
    any point in insidePoints to any point in outsidePoints.

Usage
    Example of function object specification:
    \verbatim
    leakFind
    {
        type            sets;

        writeControl    timeStep;
        interpolationScheme cell;
        setFormat       vtk;

        sets
        (
            leakFind
            {
                type    shortestPath;
                insidePoints   ((0.08 -0.020  -0.005) (-0.05 -0.020  -0.005));
                outsidePoints  ((-0.08 -0.020  -0.005)(0.05 -0.020  -0.005));
                axis    xyz;
            }
        );

        // Needs at least one field
        fields          ( p );
    }
    \endverbatim

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | shortestPath                            | yes      |
        axis     | x, y, z, xyz, distance                  | yes      |
        insidePoints  | The inside points                  | yes      |
        outsidePoints | The outside points                 | yes      |
    \endtable

SourceFiles
    shortestPathSet.C

\*---------------------------------------------------------------------------*/

#ifndef shortestPathSet_H
#define shortestPathSet_H

#include "sampledSet.H"
#include "bitSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type> class topoDistanceData;

/*---------------------------------------------------------------------------*\
                       Class shortestPathSet Declaration
\*---------------------------------------------------------------------------*/

class shortestPathSet
:
    public sampledSet
{
    // Private data

        //- Originating set of points
        const pointField insidePoints_;

        //- Destination set of points
        const pointField outsidePoints_;

        //- Local faces that were closed
        labelList leakFaces_;


    // Private Member Functions

        //- Get face with least distance along route
        static label findMinFace
        (
            const polyMesh& mesh,
            const label cellI,
            const List<topoDistanceData<label>>& allFaceInfo,
            const bitSet& isLeakPoint,
            const bool minDistance,
            const point& origin
        );

        //- Sync the leak path data
        void sync
        (
            const polyMesh& mesh,
            bitSet& isLeakFace,
            bitSet& isLeakPoint,
            const label celli,
            point& origin,
            bool& findMinDistance
        ) const;

        //- Calculate (topological) distance to cellI
        void calculateDistance
        (
            const label iter,
            const polyMesh& mesh,
            const label cellI,

            List<topoDistanceData<label>>& allFaceInfo,
            List<topoDistanceData<label>>& allCellInfo
        ) const;

        //- Checks if face uses a leak point
        bool touchesWall
        (
            const polyMesh& mesh,
            const label facei,

            bitSet& isLeakFace,
            const bitSet& isLeakPoint
        ) const;

        //- Removes open-edged faces from set. Ignores edges with both points
        //  in isBlockedPoint. Returns global number of faces removed from set.
        label erodeFaceSet
        (
            const polyMesh& mesh,
            const bitSet& isBlockedPoint,
            bitSet& isLeakFace
        ) const;

        //- Mark faces inbetween leak cells
        bitSet pathFaces(const polyMesh& mesh, const bitSet& isLeakCell) const;

        //- Calculate path between insideCelli (-1 if not on current processor)
        //  and outsidePoint. Appends cellcentres on path to track.
        //      isLeakCell  : track passes through cell
        //      isLeakFace  : faces of leak cells that are also on boundary
        //      isLeakPoint : points  of leak faces       ,,
        // Returns true if new path found, false otherwise
        bool genSingleLeakPath
        (
            const bool markLeakPath,
            const label iter,
            const polyMesh& mesh,
            const bitSet& isBlockedFace,
            const point& insidePoint,
            const label insideCelli,
            const point& outsidePoint,
            const label outsideCelli,

            // Generated sampling points
            const scalar trackLength,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist,

            // State of current leak paths
            bitSet& isLeakCell,
            bitSet& isLeakFace,
            bitSet& isLeakPoint,

            // Work storage
            List<topoDistanceData<label>>& allFaceInfo,
            List<topoDistanceData<label>>& allCellInfo
        ) const;

        //- Calculate path between insideCelli (-1 if not on current processor)
        //  and outsidePoint. Appends cellcentres on path to track.
        //      isLeakCell  : track passes through cell
        //      isLeakFace  : faces of leak cells that are also on boundary
        //      isLeakPoint : points  of leak faces       ,,
        void genSamples
        (
            const bool markLeakPath,
            const label maxIter,
            const polyMesh& mesh,
            const bitSet& isBlockedFace,
            const point& insidePoint,
            const label insideCelli,
            const point& outsidePoint,

            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist,
            bitSet& isLeakCell,
            bitSet& isLeakFace,
            bitSet& isLeakPoint
        ) const;

        //- Generate whole path. With markLeakPath it will block all faces
        //  along the whole path so will maximise the chances of finding
        //  multiple gaps. With markLeakPath=false it will only block any
        //  faces connected to a boundary. This makes for the nicest
        //  hole-filling.
        void genSamples
        (
            const bool markLeakPath,    // mark all cells along path
            const label maxIter,
            const polyMesh& mesh,
            const labelUList& wallPatches,
            const bitSet& blockedFace
        );


public:

    //- Runtime type information
    TypeName("shortestPath");


    // Constructors

        //- Construct from components. blockedFace is an optional specification
        //  of face that behave as if a wall
        shortestPathSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const bool markLeakPath,
            const label maxIter,
            const labelUList& wallPatches,
            const pointField& insidePoints,
            const pointField& outsidePoints,
            const boolList& blockedFace
        );

        //- Construct from dictionary
        shortestPathSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~shortestPathSet() = default;


    // Member Functions

        //- Set of mesh faces needed to close the gap. Will close the gap but
        //  might be bigger than required
        const labelList& leakFaces() const
        {
            return leakFaces_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
