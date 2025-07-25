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
    Foam::one

Description
    A class representing the concept of 1 (one) that can be used to avoid
    manipulating objects known to be \em one at compile-time.
    It is also used for tagged dispatch.

SourceFiles
    oneI.H

SeeAlso
    Foam::zero

\*---------------------------------------------------------------------------*/

#ifndef Foam_one_H
#define Foam_one_H

#include "labelFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class Istream;
class Ostream;

/*---------------------------------------------------------------------------*\
                             Class one Declaration
\*---------------------------------------------------------------------------*/

class one
{
public:

    //- The value type
    typedef one value_type;

    //- Default construct
    constexpr one() noexcept {}

    //- Construct from Istream consumes no content.
    explicit constexpr one(Istream&) noexcept {}


    //- Return 1 for label
    constexpr operator label() const noexcept { return 1; }

    //- Return 1 for float
    constexpr operator float() const noexcept { return 1; }

    //- Return 1 for double
    constexpr operator double() const noexcept { return 1; }

    //- Component-wise or element-wise access returns one
    one operator[](const label) const noexcept { return one{}; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Global one (1)
static constexpr const one One;

// IOstream Operators

//- Read from Istream consumes no content.
inline constexpr Istream& operator>>(Istream& is, one&) noexcept
{
    return is;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Operators, Functions
#include "oneI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
