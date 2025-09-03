/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Foam::objectHit

Description
    This class describes a combination of target object index and success flag.
    Behaves somewhat like std::optional

\*---------------------------------------------------------------------------*/

#ifndef Foam_objectHit_H
#define Foam_objectHit_H

#include "bool.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectHit;
inline Ostream& operator<<(Ostream& os, const objectHit& obj);

/*---------------------------------------------------------------------------*\
                           Class objectHit Declaration
\*---------------------------------------------------------------------------*/

class objectHit
{
    // Private Data

        //- Hit success
        bool hit_;

        //- The hit object index
        label index_;


public:

    // Constructors

        //- Default construct. Nothing hit and object index = -1
        constexpr objectHit() noexcept
        :
            hit_(false),
            index_(-1)
        {}

        //- Construct from components
        objectHit(const bool success, const label index) noexcept
        :
            hit_(success),
            index_(index)
        {}

        //- Construct from Istream
        explicit objectHit(Istream& is)
        :
            hit_(readBool(is)),
            index_(readLabel(is))
        {}


    // Member Functions

    // Access

        //- Is there a hit?
        bool hit() const noexcept
        {
            return hit_;
        }

        //- Return the hit object index
        label index() const noexcept
        {
            return index_;
        }

        //- Identical to index()
        label hitObject() const noexcept
        {
            return index_;
        }


    // Edit

        //- Reset to default construct state (false, -1)
        void reset() noexcept
        {
            hit_ = false;
            index_ = -1;
        }

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

        //- Set the hit object index
        void setIndex(const label index) noexcept
        {
            index_ = index;
        }


    // Ostream Operator

        inline friend Ostream& operator<<(Ostream& os, const objectHit& obj)
        {
             return os << obj.hit() << token::SPACE << obj.index();
        }
};


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

inline bool operator==(const objectHit& a, const objectHit& b)
{
    return a.hit() == b.hit() && a.index() == b.index();
}


inline bool operator!=(const objectHit& a, const objectHit& b)
{
    return !(a == b);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
