/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::patchEdgeSet

Description
    Like Foam::uniformSet but samples patch edges

Usage
    Example of function object specification:
    \verbatim
    sets
    {
        type                sets;
        libs                (sampling);
        writeControl        timeStep;
        writeInterval       1;

        fields              (p U);
        interpolationScheme cellPoint;
        setFormat           vtk;

        sets
        (
            // Intersections of patches with plane
            patchEdge
            {
                type        patchEdge;
                axis        x;

                // List of patches to sample (into single file). Note: supports
                // wildcards.
                patches     (movingWall);

                // Surface type
                surfaceType searchablePlane;
                // Additional info for surface
                planeType   pointAndNormal;
                pointAndNormalDict
                {
                    point   (1.5 1.5 1.5);
                    normal  (0.1 0 1);
                }

                // Sort point according to distance to origin
                origin      (0 1 0);
            }

            // Intersections of patches with stl
            sphere.stl
            {
                type        patchEdge;
                axis        x;

                patches     (movingWall);

                surfaceType triSurfaceMesh;
                // Sort point according to distance to origin
                origin      (0 1 0);
            }
        );
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property  | Description                         | Required | Default
        type      | patchEdge                           | yes      |
        axis      | x, y, z, xyz, distance              | yes      |
        patches   | List of patch names or regexs       | yes      |
        surfaceType | Definition of the surface         | yes      |
        origin    | reference location                  | yes      |
    \endtable

Note
    - The ordering of the points is according to the distance to the specified
    origin.
    - Collocated points get merged if closer than SMALL (1e-15 in double prec).
    (hardcoded). Collocated points can occur if patch points are exactly
    on the sampling plane. Note that (in parallel) any point merging gets
    done locally before sending to the master.

SourceFiles
    patchEdgeSet.C

\*---------------------------------------------------------------------------*/

#ifndef patchEdgeSet_H
#define patchEdgeSet_H

#include "sampledSet.H"
#include "DynamicList.H"
#include "HashSet.H"
#include "plane.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurface;

/*---------------------------------------------------------------------------*\
                        Class patchEdgeSet Declaration
\*---------------------------------------------------------------------------*/

class patchEdgeSet
:
    public sampledSet
{
    // Private data

        //- The surface
        const autoPtr<searchableSurface> surfPtr_;

        //- Reference point
        const point origin_;

        //- Patches to sample
        const labelHashSet patchSet_;


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
    TypeName("patchEdge");


    // Constructors

        //- Construct from dictionary
        patchEdgeSet
        (
            const word& name,
            const polyMesh& mesh,
            const meshSearch& searchEngine,
            const dictionary& dict
        );


    //- Destructor
    virtual ~patchEdgeSet() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
