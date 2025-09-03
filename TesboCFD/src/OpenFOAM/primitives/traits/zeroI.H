/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "zero.H"
#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline constexpr zero operator+(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr const Type& operator+(const Type& val, const zero&) noexcept
{
    return val;
}

template<class Type>
inline constexpr const Type& operator+(const zero&, const Type& val) noexcept
{
    return val;
}

inline constexpr zero operator-(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr const Type& operator-(const Type& val, const zero&) noexcept
{
    return val;
}

template<class Type>
inline constexpr Type operator-(const zero&, const Type& val)
{
    return -val;
}

inline constexpr zero operator*(const zero&, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator*(const Type& val, const zero&) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator*(const zero&, const Type& val) noexcept
{
    return Zero;
}

template<class Type>
inline constexpr zero operator/(const zero&, const Type& val) noexcept
{
    return Zero;
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

inline zero min(const zero&, const zero&)
{
    return Zero;
}

template<class Type>
inline Type min(const zero&, const Type& t)
{
    return min(scalar(0), t);
}

template<class Type>
inline Type min(const Type& t, const zero&)
{
    return min(t, scalar(0));
}

inline zero max(const zero&, const zero&)
{
    return Zero;
}

template<class Type>
inline Type max(const zero&, const Type& t)
{
    return max(scalar(0), t);
}

template<class Type>
inline Type max(const Type& t, const zero&)
{
    return max(t, scalar(0));
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
