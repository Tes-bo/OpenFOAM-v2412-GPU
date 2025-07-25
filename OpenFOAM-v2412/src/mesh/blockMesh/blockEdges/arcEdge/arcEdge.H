/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::blockEdges::arcEdge

Description
    A blockEdge defined as an arc of a circle.

    The arc is normally defined by its endpoints and a point on
    its circumference, typically a midpoint. For example,
    \verbatim
    points ((1 0 0) (0 1 0));

    arc 0 1 (0.707107 0.707107 0);
    \endverbatim
    The arc can enclose an angle greater than 0 and less than 360 degrees.

    The arc will frequently enclose an angle less than 180 degrees.
    For the case, it is possible to define the arc by its endpoints and its
    centre (origin) point. For example,
    \verbatim
    arc 0 1 origin (0 0 0);
    \endverbatim
    When defined in the way, any discrepancy in the arc radius for the
    endpoints is resolved by adjusting the origin to ensure that the average
    radius is satisfied.

    It is also possible to define a \em flatness factor as a multiplier of
    the calculated radius. For example,
    \verbatim
    arc 0 1 origin 1.1 (0 0 0);
    \endverbatim

SourceFiles
    arcEdge.C

\*---------------------------------------------------------------------------*/

#ifndef blockEdges_arcEdge_H
#define blockEdges_arcEdge_H

#include "blockEdge.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace blockEdges
{

/*---------------------------------------------------------------------------*\
                           Class arcEdge Declaration
\*---------------------------------------------------------------------------*/

class arcEdge
:
    public blockEdge
{
    // Private Data

        //- The arc radius
        scalar radius_;

        //- The arc angle (radians)
        scalar angle_;

        //- The local cylindrical coordinate system
        coordSystem::cylindrical cs_;


    // Private Member Functions

        //- Calculate angle, radius, cylindrical coordinate system
        //- from end points and the given point on the circumference
        void calcFromMidPoint
        (
            const point& p1,        //!< Start point
            const point& p3,        //!< End point
            const point& p2         //!< Point on circumference
        );

        //- Calculate angle, radius, cylindrical coordinate system
        //- from end points and the given origin.
        //  Optionally adjust centre to accommodate deviations in the
        //  effective radius to the end-points
        void calcFromCentre
        (
            const point& p1,        //!< Start point
            const point& p3,        //!< End point
            const point& centre,    //!< Centre
            bool adjustCentre = false,  //!< Adjust centre
            scalar rMultiplier = 1      //!< Adjust radius by this factor
        );

        //- No copy construct
        arcEdge(const arcEdge&) = delete;

        //- No copy assignment
        void operator=(const arcEdge&) = delete;


public:

    //- Runtime type information
    TypeName("arc");


    // Constructors

        //- Construct from components, given the origin of the circle
        arcEdge
        (
            const pointField& points,   //!< Referenced point field
            const point& origin,        //!< The origin of the circle
            const edge& fromTo          //!< Start/end in point field
        );

        //- Construct from components, using a point on the circumference
        arcEdge
        (
            const pointField& points,   //!< Referenced point field
            const edge& fromTo,         //!< Start/end in point field
            const point& midPoint       //!< A point on the circumference
        );

        //- Construct from components, given the origin of the circle
        arcEdge
        (
            const pointField& points,   //!< Referenced point field
            const point& origin,        //!< The origin of the circle
            const label from,           //!< Start point in point field
            const label to              //!< End point in point field
        );

        //- Construct from components, using a point on the circumference
        arcEdge
        (
            const pointField& points,   //!< Referenced point field
            const label from,           //!< Start point in point field
            const label to,             //!< End point in point field
            const point& midPoint       //!< A point on the circumference
        );

        //- Construct from Istream and point field.
        //  The Istream can specify either a point on the circumference,
        //  or with a tag to specify the origin.
        arcEdge
        (
            const dictionary& dict,
            const label index,
            const searchableSurfaces&  /*unused*/,
            const pointField& points,  //!< Referenced point field
            Istream& is
        );


    //- Destructor
    virtual ~arcEdge() = default;


    // Member Functions

        //- The point corresponding to the curve parameter [0-1]
        point position(const scalar lambda) const;

        //- The length of the curve
        scalar length() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace blockEdges
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
