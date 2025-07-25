/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::polyLine

Description
    A series of straight line segments, which can also be interpreted as
    a series of control points for splines, etc.

    A future implementation could also handle a closed polyLine.

SourceFiles
    polyLine.C

\*---------------------------------------------------------------------------*/

#ifndef polyLine_H
#define polyLine_H

#include "pointField.H"
#include "scalarList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class polyLine Declaration
\*---------------------------------------------------------------------------*/

class polyLine
{
protected:

    // Protected Data

        //- The control points or ends of each segments
        pointField points_;

        //- The real (total) line length
        scalar lineLength_;

        //- The rational (0-1) cumulative parameter value for each point
        scalarList param_;


    // Protected Member Functions

        //- Precalculate the rational cumulative parameter value
        //- and the line-length
        void calcParam();

        //- Return the line segment and the local parameter [0..1]
        //- corresponding to the global lambda [0..1]
        label localParameter(scalar& lambda) const;


public:

    // Constructors

        //- Construct from components
        explicit polyLine
        (
            const pointField& points,   //!< The poly-line points
            const bool notImplementedClosed = false
        );

        //- Construct from begin, intermediate, end points
        polyLine
        (
            const point& start,             //!< The start point
            const pointField& intermediate, //!< The intermediate points
            const point& end,               //!< The end point
            const bool notImplementedClosed = false
        );


    // Static Member Functions

        //- Concatenate begin, intermediate and end points
        static tmp<pointField> concat
        (
            const point& start,             //!< The start point
            const pointField& intermediate, //!< The intermediate points
            const point& end                //!< The end point
        );


    // Member Functions

        //- Return const-access to the control-points
        const pointField& points() const noexcept;

        //- The number of line segments
        label nSegments() const noexcept;

        //- The point position corresponding to the curve parameter
        //  0 <= lambda <= 1
        point position(const scalar) const;

        //- The point position corresponding to the local parameter
        //  0 <= lambda <= 1 on the given segment
        point position(const label segment, const scalar) const;

        //- The length of the curve
        scalar length() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
