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

Class
    Foam::endian

Description
    Help with architecture-specific aspects, primarily used
    as a namespace, but defined as a class to allow pTraits specialization.

SourceFiles
    foamEndianI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_foamEndian_H
#define Foam_foamEndian_H

#include "foamEndianFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class endian Declaration
\*---------------------------------------------------------------------------*/

class endian
{
public:

    // Public Methods

        // FUTURE?
        // #if defined (WM_BIG_ENDIAN)
        // //- The endian name as a string
        // static constexpr const char* name() noexcept { return "MSB"; }
        // #else
        // //- The endian name as a string
        // static constexpr const char* name() noexcept { return "LSB"; }
        // #endif

        //- Runtime check for big endian.
        inline static bool isBig() noexcept;

        //- Runtime check for little endian.
        inline static bool isLittle() noexcept;

        //- Byte endian swapping for 32-bits
        inline static uint32_t swap32(uint32_t);

        //- Byte endian swapping for 64-bits
        inline static uint64_t swap64(uint64_t);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamEndianI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
