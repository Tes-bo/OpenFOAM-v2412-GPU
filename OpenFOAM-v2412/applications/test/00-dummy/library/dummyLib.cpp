/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "dummyLib.H"

// We know that our options file has properly defined types here

#undef SOLVE_SIZE

#if defined WM_SP
# define PRECISION    "SP"
# define SCALAR_SIZE  (8*sizeof(float))
#elif defined(WM_SPDP)
# define PRECISION    "SPDP"
# define SCALAR_SIZE  (8*sizeof(float))
# define SOLVE_SIZE   (8*sizeof(double))
#elif defined WM_DP
# define PRECISION    "DP"
# define SCALAR_SIZE  (8*sizeof(double))
#else
# define PRECISION    "QP"
# define SCALAR_SIZE  (8*sizeof(long double))
#endif

// Test additional exported symbols
#ifdef _WIN32
    #define defineWindowsLibEntryPoint(libName)                               \
        extern "C" void lib_##libName##_entry_point() {}
#else
    #define defineWindowsLibEntryPoint(libName)  /* Nothing */
#endif


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// The 'extern C' export is independent of namespace
namespace Foam
{
    defineWindowsLibEntryPoint(dummyLib);
}

const std::string Foam::Detail::dummyLib::arch(WM_ARCH);

const std::string Foam::Detail::dummyLib::compiler(WM_COMPILER);

const std::string Foam::Detail::dummyLib::precision(PRECISION);

const int Foam::Detail::dummyLib::label_size(WM_LABEL_SIZE);

const int Foam::Detail::dummyLib::scalar_size(SCALAR_SIZE);

const int Foam::Detail::dummyLib::solveScalar_size
(
    #ifdef SOLVE_SIZE
    SOLVE_SIZE
    #else
    SCALAR_SIZE
    #endif
);

const std::string Foam::Detail::dummyLib::archComp
(
    WM_ARCH WM_COMPILER
);

const std::string Foam::Detail::dummyLib::archCompBase
(
    WM_ARCH WM_COMPILER PRECISION "Int"
  + std::to_string(WM_LABEL_SIZE)
);

const std::string Foam::Detail::dummyLib::archCompFull
(
    WM_ARCH WM_COMPILER PRECISION "Int"
  + std::to_string(WM_LABEL_SIZE)
  + WM_COMPILE_OPTION
);


// ************************************************************************* //
