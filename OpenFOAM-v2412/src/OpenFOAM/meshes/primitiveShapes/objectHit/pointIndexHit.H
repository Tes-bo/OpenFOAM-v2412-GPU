/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    Foam::PointIndexHit

Description
    This class describes the interaction of an object (often a face)
    and a point.
    It carries the info of a successful hit and (if successful),
    returns the interaction point.

    Like pointHit but carries object (eg, face, cell, edge etc.) index

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointIndexHit_H
#define Foam_pointIndexHit_H

#include "pointHit.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class PointType> class PointIndexHit;

template<class PointType>
inline Istream& operator>>(Istream& is, PointIndexHit<PointType>& pHit);

template<class PointType>
inline Ostream& operator<<(Ostream& os, const PointIndexHit<PointType>& pHit);


// Standard Types

//- A PointIndexHit with a 3D point
typedef PointIndexHit<point> pointIndexHit;


/*---------------------------------------------------------------------------*\
                        Class PointIndexHit Declaration
\*---------------------------------------------------------------------------*/

template<class PointType>
class PointIndexHit
{
    // Private Data

        //- Hit success
        bool hit_;

        //- Point of hit; invalid for misses
        PointType point_;

        //- Label of object hit
        label index_;


public:

    // Public Typedefs

        //- The point type
        typedef PointType point_type;


    // Constructors

        //- Default construct. A zero point, with no hit and index = -1
        PointIndexHit()
        :
            hit_(false),
            point_(Zero),
            index_(-1)
        {}

        //- Construct from a point, with no hit and index = -1
        explicit PointIndexHit(const point_type& p)
        :
            hit_(false),
            point_(p),
            index_(-1)
        {}

        //- Construct from hitPoint with index = -1 or specified
        explicit PointIndexHit(const PointHit<point_type>& p, label index = -1)
        :
            hit_(p.hit()),
            point_(p.point()),
            index_(index)
        {}

        //- Copy construct with different index
        PointIndexHit(const PointIndexHit<point_type>& p, label index)
        :
            hit_(p.hit()),
            point_(p.point()),
            index_(index)
        {}

        //- Construct from components
        PointIndexHit(bool success, const point_type& p, label index)
        :
            hit_(success),
            point_(p),
            index_(index)
        {}

        //- Construct from Istream
        explicit PointIndexHit(Istream& is)
        {
            is >> *this;
        }


    // Member Functions

    // Access

        //- Is there a hit?
        bool hit() const noexcept
        {
            return hit_;
        }

        //- Return the hit index
        label index() const noexcept
        {
            return index_;
        }

        //- Return point, no checks
        const point_type& point() const noexcept
        {
            return point_;
        }

        //- Access the point, no checks
        point_type& point() noexcept
        {
            return point_;
        }

        //- Return hit point. Fatal if not hit.
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

        //- Return miss point. Fatal if hit.
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

        //- The point, no checks. Same as point()
        //  \deprecated(2020-10) use point()
        const point_type& rawPoint() const noexcept { return point_; }

        //- The point, no checks. Same as point()
        //  \deprecated(2020-10) use point()
        point_type& rawPoint() noexcept { return point_; }


    // Edit

        //- Set the hit status \em on
        void setHit() noexcept
        {
            hit_ = true;
        }

        //- Set the hit status \em off
        void setMiss() noexcept
        {
            hit_ = false;
        }

        //- Set the point
        void setPoint(const point_type& p)
        {
            point_ = p;
        }

        //- Set the index
        void setIndex(const label index) noexcept
        {
            index_ = index;
        }

        //- Set the point as \em hit without changing the hit-index
        void hitPoint(const point_type& p)
        {
            point_ = p;
            hit_ = true;
        }

        //- Set the point as \em hit and set the hit-index
        void hitPoint(const point_type& p, const label index)
        {
            point_ = p;
            hit_ = true;
            index_ = index;
        }


    // Write

        //- Report hit/miss status, point and index
        void write(Ostream& os)
        {
            os  << (hit_ ? "hit:" : "miss:")
                << point_ << " index:" << index_;
        }


    // Member Operators

        //- Test for equality of all components
        bool operator==(const PointIndexHit& rhs) const
        {
            return
            (
                hit_ == rhs.hit_
             && index_ == rhs.index_
             && point_ == rhs.point_
            );
        }

        //- Test for inequality of components
        bool operator!=(const PointIndexHit& rhs) const
        {
            return !(*this == rhs);
        }


    // IO Operators

        friend Ostream& operator<<(Ostream& os, const PointIndexHit& pHit)
        {
            if (os.format() == IOstreamOption::BINARY)
            {
                os.write
                (
                    reinterpret_cast<const char*>(&pHit),
                    sizeof(PointIndexHit)
                );
            }
            else
            {
                os  << pHit.hit_ << token::SPACE
                    << pHit.point_ << token::SPACE
                    << pHit.index_;
            }

            os.check(FUNCTION_NAME);
            return os;
        }


        friend Istream& operator>>(Istream& is, PointIndexHit& pHit)
        {
            if (is.format() == IOstreamOption::BINARY)
            {
                // TODO (2019-08-06):
                // cannot properly handle mixed-precision reading
                // owing to bool and Point type.

                is.read
                (
                    reinterpret_cast<char*>(&pHit),
                    sizeof(PointIndexHit)
                );
            }
            else
            {
                is >> pHit.hit_ >> pHit.point_ >> pHit.index_;
            }

            is.check(FUNCTION_NAME);
            return is;
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for pointIndexHit
template<> struct is_contiguous<pointIndexHit> : is_contiguous<point> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
