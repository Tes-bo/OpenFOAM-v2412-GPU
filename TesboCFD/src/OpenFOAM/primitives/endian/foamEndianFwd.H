/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Description
    Define WM_BIG_ENDIAN or WM_LITTLE_ENDIAN
    and forward declare Foam::endian class.

\*---------------------------------------------------------------------------*/

#ifndef Foam_foamEndianFwd_H
#define Foam_foamEndianFwd_H

#include <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef __BYTE_ORDER__
    // Clang, Gcc, Icc, Pgi
    #if (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
        #define WM_LITTLE_ENDIAN
    #elif (__BYTE_ORDER__ == __ORDER_BIG_ENDIAN__)
        #define WM_BIG_ENDIAN
    #else
        #error "__BYTE_ORDER__ is not BIG or LITTLE endian"
    #endif
#endif

// Could downgrade to a warning, but then user always needs runtime check.
#if !defined(WM_BIG_ENDIAN) && !defined(WM_LITTLE_ENDIAN)
    #error "Cannot determine BIG or LITTLE endian."
    #error "Please add to compilation options"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class endian;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
