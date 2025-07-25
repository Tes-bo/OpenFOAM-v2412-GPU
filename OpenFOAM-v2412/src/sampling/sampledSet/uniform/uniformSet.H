/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::uniformSet

Description
    A sampler type which provides a uniform distribution of \c nPoints
    sample locations along a straight line specified between a given
    \c start and an \c end points.

Usage
    Example specification:
    \verbatim
    sets
    (
        <set>
        {
            // Mandatory entries
            type        uniform;
            axis        <options>;
            start       <vector>;
            end         <vector>;
            nPoints     <label>;

            // Optional entries
            tol         <scalar>;
        }

    );
    \endverbatim

    where the entries mean:
    \table
      Property  | Description                        | Type | Reqd  | Deflt
      type      | Type name: uniform                 | word | yes   | -
      axis      | Output type of sample locations    | word | yes   | -
      start     | Start point of sample line         | vector | yes | -
      end       | End point of sample line           | vector | yes | -
      nPoints   | Number of points between start/end | label  | yes | -
      tol       | Relative tolerance                 | scalar | no  | 1e-3
    \endtable

    Options for the \c axis entry:
    \verbatim
      x        | x-ordinate of a sample
      y        | y-ordinate of a sample
      z        | z-ordinate of a sample
      xyz      | x-y-z coordinates of a sample
      distance | Normal distance to the first point of a sample
    \endverbatim

SourceFiles
    uniformSet.C

\*---------------------------------------------------------------------------*/

#ifndef uniformSet_H
#define uniformSet_H

#include "passiveParticleCloud.H"
#include "sampledSet.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class uniformSet Declaration
\*---------------------------------------------------------------------------*/

class uniformSet
:
    public sampledSet
{
    // Private Data

        //- Start point of sample line
        point start_;

        //- End point of sample line
        point end_;

        //- Relative tolerance when comparing points
        //- relative to difference between start_ and end_
        scalar tol_;

        //- Number of sampling points
        label nPoints_;


    // Private Member Functions

        //- Calculates - starting at samplePt - the first sampling point
        //  on or after currentPt. smallDist is the tolerance used to compare
        //  positions. Returns false if end of samples reached.
        bool nextSample
        (
            const point& currentPt,
            const vector& offset,
            const scalar smallDist,
            point& samplePt,
            label& sampleI
        ) const;

        //- Samples from startTrackPt/Celli. Updates particle/samplePt/sampleI
        //  and puts
        //  samples in the DynamicLists. Returns false if end of all samples
        //  reached
        bool trackToBoundary
        (
            passiveParticleCloud& particleCloud,
            passiveParticle& singleParticle,
            point& samplePt,
            label& sampleI,
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Samples from start_ to end_. samplingSegments contains segmentNo
        //  for each sample.
        void calcSamples
        (
            DynamicList<point>& samplingPts,
            DynamicList<label>& samplingCells,
            DynamicList<label>& samplingFaces,
            DynamicList<label>& samplingSegments,
            DynamicList<scalar>& samplingCurveDist
        ) const;

        //- Uses calcSamples to obtain samples. Copies them into *this.
        void genSamples();


public:

    //- Runtime type information
    TypeName("uniform");


    // Constructors

        //- Construct from components
        uniformSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const point& start,
            const point& end,
            const label nPoints
        );

        //- Construct from dictionary
        uniformSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~uniformSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
