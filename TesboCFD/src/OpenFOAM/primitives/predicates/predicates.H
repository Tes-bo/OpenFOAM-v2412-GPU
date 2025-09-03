/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenCFD Ltd.
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

Namespace
    Foam::predicates

Description
    Constant predicate types.

SourceFiles
    predicates.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_predicates_H
#define Foam_predicates_H

#include <string>
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace predicates
{

/*---------------------------------------------------------------------------*\
                           Class always Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates that always return true, useful for templating.
struct always : std::true_type
{
    //- Unary predicate
    //  \return true
    template<class T>
    constexpr bool operator()(const T&) const noexcept
    {
        return value;
    }

    //- Binary predicate
    //  \return true
    template<class T1, class T2>
    constexpr bool operator()(const T1&, const T2&) const noexcept
    {
        return value;
    }

    //- String match
    //  \return true
    constexpr bool match(const std::string&, bool literal=false) const noexcept
    {
        return value;
    }
};


/*---------------------------------------------------------------------------*\
                            Class never Declaration
\*---------------------------------------------------------------------------*/

//- Unary and binary predicates that never return true, useful for templating.
struct never : std::false_type
{
    //- Unary predicate
    //  \return false
    template<class T>
    constexpr bool operator()(const T&) const noexcept
    {
        return value;
    }

    //- Binary predicate
    //  \return false
    template<class T1, class T2>
    constexpr bool operator()(const T1&, const T2&) const noexcept
    {
        return value;
    }

    //- String match
    //  \return false
    constexpr bool match(const std::string&, bool literal=false) const noexcept
    {
        return value;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace predicates
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
