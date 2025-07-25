/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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
    Foam::zero

Description
    A class representing the concept of 0 (zero) that can be used to avoid
    manipulating objects known to be \em zero at compile-time.
    It is also used for tagged dispatch.

SourceFiles
    zero.C
    zeroI.H

SeeAlso
    Foam::one

\*---------------------------------------------------------------------------*/

#ifndef Foam_zero_H
#define Foam_zero_H

#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                            Class zero Declaration
\*---------------------------------------------------------------------------*/

class zero
{
public:

    //- The value type
    typedef zero value_type;

    //- A static instance, when modifiable reference is required by an API
    static zero dummy;

    //- Default construct
    constexpr zero() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr zero(Istream&) noexcept {}


    //- Return false (0) for bool
    constexpr operator bool() const noexcept { return false; }

    //- Return 0 for label
    constexpr operator label() const noexcept { return 0; }

    //- Return 0 for float
    constexpr operator float() const noexcept { return 0; }

    //- Return 0 for double
    constexpr operator double() const noexcept { return 0; }

    //- Component-wise or element-wise access returns zero
    zero operator[](const label) const noexcept { return zero{}; }

    // FUTURE?: Swallow assignment (as per std::ignore)
    // template<class T>
    // constexpr const zero& operator=(const T&) const { return *this; }

    //- Unary functor returns zero
    template<class T>
    constexpr zero operator()(const T&) const noexcept
    {
        return zero{};
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global zero (0)
static constexpr const zero Zero;

// IOstream Operators

//- Read from Istream consumes no content
inline constexpr Istream& operator>>(Istream& is, zero&) noexcept
{
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators, Functions

#include "zeroI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
