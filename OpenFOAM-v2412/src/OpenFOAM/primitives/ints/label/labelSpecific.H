/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

Typedef
    Foam::label

Description
    label specific functions for which there is no uLabel equivalent

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return 1 if s is greater_equal zero, or otherwise -1
inline label sign(const label s) noexcept
{
    return (s >= 0)? 1: -1;
}


//- Return 1 if s is greater_equal zero, or otherwise 0
inline label pos0(const label s) noexcept
{
    return (s >= 0)? 1: 0;
}


//- Return 1 if s is less than zero, or otherwise 0
inline label neg(const label s) noexcept
{
    return (s < 0)? 1: 0;
}


//- Return 1 if s is less_equal zero, or otherwise 0
inline label neg0(const label s) noexcept
{
    return (s <= 0)? 1: 0;
}


//- Return the positive part of s, otherwise zero. Same as max(0, s).
inline label posPart(const label s) noexcept
{
    return (s > 0)? s: 0;
}


//- Return the negative part of s, otherwise zero. Same as min(0, s).
//  Does not change the sign
inline label negPart(const label s) noexcept
{
    return (s < 0)? s: 0;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
