/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2014 OpenFOAM Foundation
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
    Foam::patchSeedSet

Description
    Initialises points on or just off patch

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | patchSeed                               | yes      |
        axis     | x, y, z, xyz, distance                  | yes      |
        patches  | List of patch names or regexs           | yes      |
        maxPoints | Max number of points to seed           | yes      |
        points   | List of selected locations              | no       | empty
    \endtable

SourceFiles
    patchSeedSet.C

\*---------------------------------------------------------------------------*/

#ifndef patchSeedSet_H
#define patchSeedSet_H

#include "sampledSet.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "Random.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class patchSeedSet Declaration
\*---------------------------------------------------------------------------*/

class patchSeedSet
:
    public sampledSet
{
    // Private data

        //- Patches to sample
        const labelHashSet patchSet_;

        //- Maximum number of patch faces to seed (if in random subset mode)
        const label maxPoints_;

        //- Random number generator (if maxPoints < num patch faces)
        autoPtr<Random> rndGenPtr_;

        //- Patch faces to seed selected based on nearness to supplied points
        const pointField selectedLocations_;


    // Private Member Functions

        //- Determine (random) samples on patch
        void calcPatchSamples
        (
            const label nAvailable,
            const label nPatchPoints,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        );

        //- Determine samples provided by supplied points
        void calcSelectedLocations
        (
            const label nAvailable,
            const label nPatchPoints,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        );

        //- Samples all points in sampleCoords.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        );

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("patchSeed");


    // Constructors

        //- Construct from dictionary
        patchSeedSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchSeedSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
