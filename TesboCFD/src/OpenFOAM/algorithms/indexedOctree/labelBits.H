/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2022 OpenCFD Ltd.
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
    Foam::labelBits

Description
    A 29bits (or 61bits) integer label with 3bits direction
    (eg, octant) packed into single label

\*---------------------------------------------------------------------------*/

#ifndef Foam_labelBits_H
#define Foam_labelBits_H

#include "label.H"
#include "direction.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class labelBits;

/*---------------------------------------------------------------------------*\
                          Class labelBits Declaration
\*---------------------------------------------------------------------------*/

class labelBits
{
    // Private Data

        //- Integer value and (octant) direction
        label data_;

public:

    // Static Functions

        #ifdef FULLDEBUG

        //- Pack integer value and bits (octant) into a label
        inline static label pack
        (
            const label val,
            const direction bits
        )
        {
            if (bits > 7 || (((val<<3)>>3) != val))
            {
                FatalErrorInFunction
                    << "Direction " << bits << " outside range 0..7"
                    << " or value " << val << " negative or larger than "
                    << label(8*sizeof(label)-3) << " bit representation"
                    << abort(FatalError);
            }
            return (val<<3) | bits;
        }

        #else

        //- Pack integer value and bits (octant) into a label
        static constexpr label pack
        (
            const label val,
            const direction bits
        ) noexcept
        {
            return (val<<3) | bits;
        }

        #endif


    // Constructors

        //- Default construct as zero initialized
        constexpr labelBits() noexcept
        :
             data_(0)
        {}

        //- Construct from components
        labelBits(const label val, const direction bits)
        :
             data_(pack(val, bits))
        {}

        //- Read construct from Istream
        explicit labelBits(Istream& is)
        {
            is >> data_;
        }


    // Member Functions

        //- The raw data value
        label data() const noexcept { return data_; }

        //- Return the integer value
        label val() const noexcept { return (data_ >> 3); }

        //- Return the octant direction
        direction bits() const noexcept { return (data_ & 7); }

        //- Set the integer value
        void setVal(const label val)
        {
            data_ = pack(val, bits());
        }

        //- Set the octant direction
        void setBits(const direction bits)
        {
            data_ = pack(val(), bits);
        }


    // Member Operators

        inline friend
        bool operator==(const labelBits& a, const labelBits& b)
        {
            return a.data_ == b.data_;
        }

        inline friend
        bool operator!=(const labelBits& a, const labelBits& b)
        {
            return !(a == b);
        }


    // IOstream Operators

        inline friend
        Istream& operator>>(Istream& is, labelBits& rhs)
        {
            return is >> rhs.data_;
        }

        inline friend
        Ostream& operator<<(Ostream& os, const labelBits& rhs)
        {
            return os << rhs.data_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
