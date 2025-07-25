/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::arraySet

Description
    Specifies an x,y,z array of uniformly distributed sampling points.

    For a dictionary specification:
    \table
        Property | Description                             | Required | Default
        type     | array                                   | yes      |
        axis     | x, y, z, xyz, distance                  | yes      |
        pointsDensity | The sampling density as (x y z) integers | yes      |
        spanBox   | The sample box dimensions (vector)     | yes      |
    \endtable

    The dictionary should also contain an embedded (Cartesian) coordinateSystem
    specification, for example,

    \verbatim
        origin (-0.05 -0.05 -0.05);
        rotation
        {
            type    axes;
            e1 (1 0 0);
            e2 (0 1 0);
        }
    \endverbatim

SourceFiles
    arraySet.C

\*---------------------------------------------------------------------------*/

#ifndef arraySet_H
#define arraySet_H

#include "sampledSet.H"
#include "labelVector.H"
#include "DynamicList.H"
#include "cartesianCS.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class arraySet Declaration
\*---------------------------------------------------------------------------*/

class arraySet
:
    public sampledSet
{
    // Private data

        //- Local Cartesian coordinate system
        coordSystem::cartesian csys_;

        //- Point density vector
        labelVector pointsDensity_;

        //- Span box
        vector spanBox_;


    // Private Member Functions

        //- Samples all points in sampleCoords.
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
    TypeName("array");


    // Constructors

        //- Construct from components
        arraySet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const word& axis,
            const coordSystem::cartesian& csys,
            const Vector<label>& pointsDensity,
            const Vector<scalar>& spanBox
        );

        //- Construct from dictionary
        arraySet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~arraySet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
