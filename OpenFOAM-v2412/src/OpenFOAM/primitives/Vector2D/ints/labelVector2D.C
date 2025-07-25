/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "labelVector2D.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelVector2D::vsType::typeName = "labelVector2D";


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#undef  defineTraits
#define defineTraits(Type, Prefix)                                            \
                                                                              \
    template<>                                                                \
    const char* const Foam::Vector2D<Type>::vsType::componentNames[] =        \
    {                                                                         \
        "x", "y"                                                              \
    };                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::zero             \
    (                                                                         \
        Vector2D<Type>::uniform(0)                                            \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::one              \
    (                                                                         \
        Vector2D<Type>::uniform(1)                                            \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::max              \
    (                                                                         \
        Vector2D<Type>::uniform(Prefix##Max)                                  \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::min              \
    (                                                                         \
        Vector2D<Type>::uniform(-Prefix##Max)                                 \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::rootMax          \
    (                                                                         \
        Vector2D<Type>::uniform(::sqrt(double(Prefix##Max)))                  \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector2D<Type> Foam::Vector2D<Type>::vsType::rootMin          \
    (                                                                         \
        Vector2D<Type>::uniform(-::sqrt(double(Prefix##Max)))                 \
    );


defineTraits(Foam::label, label);

#undef defineTraits

// ************************************************************************* //
