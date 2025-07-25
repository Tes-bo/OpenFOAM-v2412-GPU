/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 Ivor Clifford/Paul Scherrer Institut
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
    Foam::extrudeModels::polyline

Description
    Extrudes by transforming points along a polyline provided as a
    series of points and edge segments. Supports all blockMesh edge
    types, e.g. line, arc, spline. The surface points are rotated to
    follow the path.

    \table
        Property  | Description                             | Required | Default
        vertices  | List of vertices                        | yes      |
        edges     | List of blockEdge segments              | yes      |
        toleranceCheck | Relative tolerance for polyline checks | no | SMALL
    \endtable

\*---------------------------------------------------------------------------*/

#ifndef polyline_H
#define polyline_H

#include "extrudeModel.H"
#include "blockEdgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace extrudeModels
{

/*---------------------------------------------------------------------------*\
                    Class polyline Declaration
\*---------------------------------------------------------------------------*/

class polyline
:
    public extrudeModel
{
    // Private Data

        //- Dummy object needed to use blockEdge
        searchableSurfaces geometry_;

        //- List of points
        pointField vertices_;

        //- List of line segments
        blockEdgeList segments_;

        //- Relative length along all segments for interplation (0 <= x <= 1)
        scalarField x_;

        //- Relative position on segments for interpolation
        //- (0 <= y <= segments_.size())
        scalarField y_;

        //- Position vector at start of polyline
        vector p0_;

        //- Direction vector at start of polyline
        vector n0_;

        //- Relative tolerance for checking alignment of polyline and surface
        scalar relTol_;

        //- Small delta for numerical differencing
        const scalar DELTA = 1e-6;

public:

    //- Runtime type information
    TypeName("polyline");


    // Constructors

        //- Construct from dictionary
        explicit polyline(const dictionary& dict);


    //- Destructor
    virtual ~polyline() = default;


    // Member Operators

        point operator()
        (
            const point& surfacePoint,
            const vector& surfaceNormal,
            const label layer
        ) const;

        //- The point and direction vector corresponding to the polyline
        //- parameter [0-1]
        void positionAndDirection
        (
            const scalar lambda,
            vector& p,
            vector& n
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace extrudeModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
