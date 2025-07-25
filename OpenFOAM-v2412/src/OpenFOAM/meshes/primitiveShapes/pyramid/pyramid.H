/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::pyramid

Description
    A geometric pyramid primitive with a base of 'n' sides:
    i.e. a parametric pyramid. A pyramid is constructed from
    a base polygon and an apex point.

SourceFiles
    pyramidI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_pyramid_H
#define Foam_pyramid_H

#include "point.H"
#include "face.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Point, class PointRef, class PolygonRef> class pyramid;

template<class Point, class PointRef, class PolygonRef>
inline Istream& operator>>
(
    Istream& is,
    pyramid<Point, PointRef, PolygonRef>& p
);

template<class Point, class PointRef, class PolygonRef>
inline Ostream& operator<<
(
    Ostream& os,
    const pyramid<Point, PointRef, PolygonRef>& p
);


// Common Typedefs

//- A pyramid using referred point and face
typedef pyramid<point, const point&, const face&> pyramidPointFaceRef;


/*---------------------------------------------------------------------------*\
                           Class pyramid Declaration
\*---------------------------------------------------------------------------*/

template<class Point, class PointRef, class PolygonRef>
class pyramid
{
    // Private Data

        PolygonRef base_;

        PointRef apex_;


public:

    // Public Typedefs

        //- The point type
        typedef Point point_type;


    // Constructors

        //- Construct from base polygon and apex point
        inline pyramid(PolygonRef base, const Point& apex);

        //- Construct from Istream
        inline explicit pyramid(Istream& is);


    // Member Functions

    // Access

        //- The apex point
        const Point& apex() const { return apex_; }

        //- The base polygon
        PolygonRef base() const { return base_; }


    // Properties

        //- Return centre (centroid)
        inline Point centre(const UList<point>& points) const;

        //- Return height vector
        inline vector height(const UList<point>& points) const;

        //- Return scalar magnitude - returns volume of pyramid
        inline scalar mag(const UList<point>& points) const;


    // IOstream Operators

        friend Istream& operator>> <Point, PointRef, PolygonRef>
        (
            Istream& is,
            pyramid& p
        );

        friend Ostream& operator<< <Point, PointRef, PolygonRef>
        (
            Ostream& os,
            const pyramid& p
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pyramidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
