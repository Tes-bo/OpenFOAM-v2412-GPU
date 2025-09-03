/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::PointHit

Description
    Describes the interaction of a object and a (templated) point.
    It carries the info of a successful hit and (if successful)
    the interaction point.

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointHit_H
#define Foam_pointHit_H

#include "bool.H"
#include "point.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class PointType> class PointHit;


// Standard Types

//- A PointHit with a 3D point
typedef PointHit<point> pointHit;


/*---------------------------------------------------------------------------*\
                           Class PointHit Declaration
\*---------------------------------------------------------------------------*/

template<class PointType>
class PointHit
{
    // Private Data

        //- Point of hit; for miss holds best estimate outside the object
        PointType point_;

        //- Distance to hit point
        scalar distance_;

        //- Hit success
        bool hit_;

        //- Eligible miss
        bool eligibleMiss_;


public:

    // Public Typedefs

        //- The point type
        typedef PointType point_type;


    // Constructors

        //- Default construct. A zero point with a large distance,
        //- no hit, no eligible misses
        PointHit()
        :
            point_(Zero),
            distance_(GREAT),
            hit_(false),
            eligibleMiss_(false)
        {}

        //- Construct from point with a large distance,
        //- no hit, no eligible misses
        explicit PointHit(const point_type& p)
        :
            point_(p),
            distance_(GREAT),
            hit_(false),
            eligibleMiss_(false)
        {}

        //- Construct from components
        PointHit
        (
            bool hit,
            const point_type& p,
            scalar dist,
            bool eligibleMiss = false
        )
        :
            point_(p),
            distance_(dist),
            hit_(hit),
            eligibleMiss_(eligibleMiss)
        {}


    // Member Functions

    // Access

        //- Is there a hit
        bool hit() const noexcept
        {
            return hit_;
        }

        //- Is this an eligible miss
        bool eligibleMiss() const noexcept
        {
            return eligibleMiss_;
        }

        //- Return the point, no checks
        const point_type& point() const noexcept
        {
            return point_;
        }

        //- Return distance to hit
        scalar distance() const noexcept
        {
            return distance_;
        }

        //- Return the hit point. Fatal if not hit.
        const point_type& hitPoint() const
        {
            if (!hit_)
            {
                FatalErrorInFunction
                    << "Requested a hit point, but it was not hit"
                    << abort(FatalError);
            }

            return point_;
        }

        //- Return the miss point. Fatal if hit.
        const point_type& missPoint() const
        {
            if (hit_)
            {
                FatalErrorInFunction
                    << "Requested a miss point, but it was hit"
                    << abort(FatalError);
            }

            return point_;
        }

        //- The point, no checks
        //  \deprecated(2020-10) use point()
        const point_type& rawPoint() const noexcept { return point_; }


    // Edit

        //- Set the hit status \em on
        void setHit() noexcept
        {
            hit_ = true;
            eligibleMiss_ = false;
        }

        //- Set the hit status \em off and set the eligible miss status
        void setMiss(const bool eligible) noexcept
        {
            hit_ = false;
            eligibleMiss_ = eligible;
        }

        //- Set the point
        void setPoint(const point_type& p)
        {
            point_ = p;
        }

        //- Set the distance
        void setDistance(const scalar d) noexcept
        {
            distance_ = d;
        }

        //- Set the point as \em hit.
        void hitPoint(const point_type& p)
        {
            point_ = p;
            hit_ = true;
            eligibleMiss_ = false;
        }


    // Member Operators

        //- Distance comparision operator, for sorting
        template<class AnyPointType>
        bool operator<(const PointHit<AnyPointType>& rhs) const noexcept
        {
            return distance_ < rhs.distance_;
        }
};


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

// Ostream Operator

template<class PointType>
inline Ostream& operator<<(Ostream& os, const PointHit<PointType>& pHit)
{
    os  << pHit.hit() << token::SPACE
        << pHit.point() << token::SPACE
        << pHit.distance() << token::SPACE
        << pHit.eligibleMiss();

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
