/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::refCount

Description
    Reference counter for various OpenFOAM components.

See also
    Foam::tmp
    Foam::token::compound

\*---------------------------------------------------------------------------*/

#ifndef Foam_refCount_H
#define Foam_refCount_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class refCount Declaration
\*---------------------------------------------------------------------------*/

class refCount
{
    // Private Data
    int count_;


public:

    //- A non-counting (dummy) refCount
    struct zero {};


    // Constructors

        //- Default construct, initializing count to 0
        constexpr refCount() noexcept
        :
            count_(0)
        {}


    // Member Functions

        //- Return the current reference count
        int use_count() const noexcept
        {
            return count_;
        }

        //- Return true if the reference count is zero
        bool unique() const noexcept
        {
            return !count_;
        }


    // Member Operators

        //- Increment the reference count
        void operator++() noexcept
        {
            ++count_;
        }

        //- Increment the reference count
        void operator++(int) noexcept
        {
            ++count_;
        }

        //- Decrement the reference count
        void operator--() noexcept
        {
            --count_;
        }

        //- Decrement the reference count
        void operator--(int) noexcept
        {
            --count_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
