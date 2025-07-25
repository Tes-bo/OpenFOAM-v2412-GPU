/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "vector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if defined(WM_DP)

template<>
const char* const Foam::Vector<float>::vsType::typeName = "floatVector";

template<>
const char* const Foam::Vector<double>::vsType::typeName = "vector";

#else

// WM_SP, WM_SPDP
template<>
const char* const Foam::Vector<float>::vsType::typeName = "vector";

template<>
const char* const Foam::Vector<double>::vsType::typeName = "doubleVector";

// or (TDB):
//
// #if defined(WM_SPDP)
// template<>
// const char* const Foam::Vector<double>::vsType::typeName = "solveVector";
// #else
// template<>
// const char* const Foam::Vector<double>::vsType::typeName = "doubleVector";
// #endif

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  defineTraits
#define defineTraits(Type, Prefix)                                            \
                                                                              \
    template<>                                                                \
    const char* const Foam::Vector<Type>::vsType::componentNames[] =          \
    {                                                                         \
        "x", "y", "z"                                                         \
    };                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::zero                 \
    (                                                                         \
        Vector<Type>::uniform(0)                                              \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::one                  \
    (                                                                         \
        Vector<Type>::uniform(1)                                              \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::max                  \
    (                                                                         \
        Vector<Type>::uniform(Prefix##VGREAT)                                 \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::min                  \
    (                                                                         \
        Vector<Type>::uniform(-Prefix##VGREAT)                                \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::rootMax              \
    (                                                                         \
        Vector<Type>::uniform(Prefix##ROOTVGREAT)                             \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::rootMin              \
    (                                                                         \
        Vector<Type>::uniform(-Prefix##ROOTVGREAT)                            \
    );


defineTraits(float, floatScalar);
defineTraits(double, doubleScalar);

#undef defineTraits

// ************************************************************************* //
