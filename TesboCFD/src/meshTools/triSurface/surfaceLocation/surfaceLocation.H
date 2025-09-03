/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::surfaceLocation

Description
    Contains information about location on a triSurface

    Access to data:
      - pointIndexHit provides
        - location
        - bool: hit/miss
        - index (of triangle/point/edge)
      - elementType() provides
        - what index above relates to. In triangle::proxType
      - triangle() provides
        - last known triangle

SourceFiles
    surfaceLocation.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_surfaceLocation_H
#define Foam_surfaceLocation_H

#include "pointIndexHit.H"
#include "triangle.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class surfaceLocation;
class triSurface;

Istream& operator>>(Istream&, surfaceLocation&);
Ostream& operator<<(Ostream&, const surfaceLocation&);

Ostream& operator<<(Ostream&, const InfoProxy<surfaceLocation>&);


/*---------------------------------------------------------------------------*\
                       Class surfaceLocation Declaration
\*---------------------------------------------------------------------------*/

class surfaceLocation
:
    public pointIndexHit
{
    // Private Data

        triPointRef::proxType elementType_;

        label triangle_;


public:

    // Constructors

        //- Default construct
        surfaceLocation()
        :
            pointIndexHit(),
            elementType_(triPointRef::NONE),
            triangle_(-1)
        {}

        //- Construct from components
        surfaceLocation
        (
            const pointIndexHit& pHit,
            const triPointRef::proxType elementType,
            const label triangle
        )
        :
            pointIndexHit(pHit),
            elementType_(elementType),
            triangle_(triangle)
        {}

        //- Construct from Istream
        explicit surfaceLocation(Istream& is)
        :
            pointIndexHit(is),
            elementType_(triPointRef::proxType(readLabel(is))),
            triangle_(readLabel(is))
        {}


    // Member Functions

        triPointRef::proxType& elementType() noexcept
        {
            return elementType_;
        }

        triPointRef::proxType elementType() const noexcept
        {
            return elementType_;
        }

        label& triangle() noexcept
        {
            return triangle_;
        }

        label triangle() const noexcept
        {
            return triangle_;
        }

        //- Normal. Approximate for points.
        vector normal(const triSurface& s) const;

        //- Return info proxy,
        //- to print information to a stream
        InfoProxy<surfaceLocation> info() const noexcept { return *this; }

        //- Write info about selected face index to a stream
        void write(Ostream& os, const triSurface& s) const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, surfaceLocation& sl);

        friend Ostream& operator<<(Ostream& os, const surfaceLocation& sl);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
