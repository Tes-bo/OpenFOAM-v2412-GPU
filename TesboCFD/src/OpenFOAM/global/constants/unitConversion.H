/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2017-2019 OpenCFD Ltd.
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

InNamespace
    Foam

Description
    Unit conversion functions

\*---------------------------------------------------------------------------*/

#ifndef unitConversion_H
#define unitConversion_H

#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Conversion from degrees to radians
inline constexpr scalar degToRad(const scalar deg) noexcept
{
    return (deg*M_PI/180.0);
}

//- Conversion from radians to degrees
inline constexpr scalar radToDeg(const scalar rad) noexcept
{
    return (rad*180.0/M_PI);
}

//- Multiplication factor for degrees to radians conversion
inline constexpr scalar degToRad() noexcept
{
    return (M_PI/180.0);
}

//- Multiplication factor for radians to degrees conversion
inline constexpr scalar radToDeg() noexcept
{
    return (180.0/M_PI);
}


//- Conversion from revolutions/minute to radians/sec
inline constexpr scalar rpmToRads(const scalar rpm) noexcept
{
    return (rpm*M_PI/30.0);
}

//- Conversion from radians/sec to revolutions/minute
inline constexpr scalar radsToRpm(const scalar rads) noexcept
{
    return (rads*30.0/M_PI);
}

//- Multiplication factor for revolutions/minute to radians/sec
inline constexpr scalar rpmToRads() noexcept
{
    return (M_PI/30.0);
}

//- Multiplication factor for radians/sec to revolutions/minute
inline constexpr scalar radsToRpm() noexcept
{
    return (30.0/M_PI);
}


//- Conversion from atm to Pa
inline constexpr scalar atmToPa(const scalar atm) noexcept
{
    return (atm*101325.0);
}

//- Conversion from bar to Pa
inline constexpr scalar barToPa(const scalar bar) noexcept
{
    return (bar*100000.0);
}


//- Conversion from Pa to atm
inline constexpr scalar paToAtm(const scalar pa) noexcept
{
    return (pa/101325.0);
}

//- Conversion from Pa to bar
inline constexpr scalar paToBar(const scalar pa) noexcept
{
    return (pa/100000.0);
}


//- User literal for degrees to radians conversion (integers)
inline constexpr scalar operator "" _deg(unsigned long long int deg) noexcept
{
    return (deg*M_PI/180.0);
}

//- User literal for degrees to radians conversion (floats)
inline constexpr scalar operator "" _deg(long double deg) noexcept
{
    return (deg*M_PI/180.0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
