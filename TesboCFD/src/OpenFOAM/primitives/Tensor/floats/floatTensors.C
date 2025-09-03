/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "tensor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#if defined(WM_DP)

template<>
const char* const Foam::Tensor<float>::vsType::typeName = "floatTensor";

template<>
const char* const Foam::Tensor<double>::vsType::typeName = "tensor";

#else
// WM_SP, WM_SPDP

template<>
const char* const Foam::Tensor<float>::vsType::typeName = "tensor";

template<>
const char* const Foam::Tensor<double>::vsType::typeName = "doubleTensor";

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  defineTraits
#define defineTraits(Type, Prefix)                                            \
                                                                              \
    template<>                                                                \
    const char* const Foam::Tensor<Type>::vsType::componentNames[] =          \
    {                                                                         \
        "xx", "xy", "xz",                                                     \
        "yx", "yy", "yz",                                                     \
        "zx", "zy", "zz"                                                      \
    };                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::zero                 \
    (                                                                         \
        Tensor<Type>::uniform(0)                                              \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::one                  \
    (                                                                         \
        Tensor<Type>::uniform(1)                                              \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::max                  \
    (                                                                         \
        Tensor<Type>::uniform(Prefix##VGREAT)                                 \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::min                  \
    (                                                                         \
        Tensor<Type>::uniform(-Prefix##VGREAT)                                \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::rootMax              \
    (                                                                         \
        Tensor<Type>::uniform(Prefix##ROOTVGREAT)                             \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::vsType::rootMin              \
    (                                                                         \
        Tensor<Type>::uniform(-Prefix##ROOTVGREAT)                            \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Tensor<Type> Foam::Tensor<Type>::I                            \
    (                                                                         \
        1, 0, 0,                                                              \
        0, 1, 0,                                                              \
        0, 0, 1                                                               \
    );                                                                        \


defineTraits(float, floatScalar);
defineTraits(double, doubleScalar);

#undef defineTraits

// ************************************************************************* //
