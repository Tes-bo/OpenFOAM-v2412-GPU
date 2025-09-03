/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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
    Foam::volumeType

Description
    An enumeration wrapper for classification of a location as being
    inside/outside of a volume.

SourceFiles
    volumeType.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_volumeType_H
#define Foam_volumeType_H

#include "contiguous.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class volumeType;
Istream& operator>>(Istream& is, volumeType& vt);
Ostream& operator<<(Ostream& os, const volumeType& vt);


/*---------------------------------------------------------------------------*\
                         Class volumeType Declaration
\*---------------------------------------------------------------------------*/

class volumeType
{
public:

    //- Volume classification types
    enum type : char
    {
        UNKNOWN = 0,    //!< Unknown state
        INSIDE  = 0x1,  //!< A location inside the volume
        OUTSIDE = 0x2,  //!< A location outside the volume
        MIXED   = 0x3   //!< A location that is partly inside and outside
    };

    // Static data

    //- Names for the classification enumeration
    static const Enum<volumeType::type> names;


private:

    // Private data

        //- Volume type
        type t_;


public:

    // Constructors

        //- Default construct as \c UNKNOWN state
        constexpr volumeType() noexcept
        :
            t_(UNKNOWN)
        {}

        //- Implicit construct from enumeration
        volumeType(type t) noexcept
        :
            t_(t)
        {}

        //- Construct as getOrDefault by name from dictionary
        volumeType(const word& key, const dictionary& dict, const type deflt);

        //- Construct from integer
        explicit volumeType(const int t)
        :
            t_(static_cast<volumeType::type>(t & 0x3))
        {}


    // Member Functions

        //- Implicit cast to the enumeration
        operator type() const noexcept
        {
            return t_;
        }

        //- The string representation of the volume type enumeration
        const word& str() const;


    // IOstream Operators

        friend Istream& operator>>(Istream& is, volumeType& vt);
        friend Ostream& operator<<(Ostream& os, const volumeType& vt);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for volumeType
template<> struct is_contiguous<volumeType> : std::true_type {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
