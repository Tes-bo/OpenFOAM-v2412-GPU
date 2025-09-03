/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "labelVector.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* const Foam::labelVector::vsType::typeName = "labelVector";

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
        Vector<Type>::uniform(Prefix##Max)                                    \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::min                  \
    (                                                                         \
        Vector<Type>::uniform(-Prefix##Max)                                   \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::rootMax              \
    (                                                                         \
        Vector<Type>::uniform(::sqrt(double(Prefix##Max)))                    \
    );                                                                        \
                                                                              \
    template<>                                                                \
    const Foam::Vector<Type> Foam::Vector<Type>::vsType::rootMin              \
    (                                                                         \
        Vector<Type>::uniform(-::sqrt(double(Prefix##Max)))                   \
    );


defineTraits(Foam::label, label);

#undef defineTraits

// ************************************************************************* //
