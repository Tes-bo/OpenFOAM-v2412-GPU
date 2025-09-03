/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include <cstdio>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Could also consider generalizing with C++11 variadic templates

template<class PrimitiveType>
std::string::size_type Foam::string::string_printf
(
    std::string& output,
    const char* fmt,
    const PrimitiveType& val
)
{
    // Use snprintf with zero to establish the size (without '\0') required
    int n = ::snprintf(nullptr, 0, fmt, val);
    if (n > 0)
    {
        output.resize(n+1);
        char* buf = &(output[0]);

        // Print directly into buffer, no stripping desired
        n = ::snprintf(buf, n+1, fmt, val);
        output.resize(n);
    }
    else
    {
        output.clear();
    }

    return output.size();
}


template<class PrimitiveType>
std::string::size_type Foam::string::string_printf
(
    std::string& output,
    const std::string& fmt,
    const PrimitiveType& val
)
{
    return string_printf(output, fmt.c_str(), val);
}


// ************************************************************************* //
