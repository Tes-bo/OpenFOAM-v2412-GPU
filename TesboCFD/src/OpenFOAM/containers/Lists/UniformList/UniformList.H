/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::UniformList

Description
    A single value that is represented as a list with an
    operator[] to access the value.
    This can be useful for templated operations expecting a list accessor.

Note
    The list currently has no sizing associated with it.

\*---------------------------------------------------------------------------*/

#ifndef Foam_UniformList_H
#define Foam_UniformList_H

#include "labelFwd.H"
#include <utility>  // For std::move

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class UniformList Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class UniformList
{
    // Private Data

        //- The value to be returned.
        Type value_;

public:

    // Constructors

        //- Construct from given value
        explicit UniformList(const Type& val) noexcept
        :
            value_(val)
        {}

        //- Move construct from given value
        explicit UniformList(Type&& val) noexcept
        :
            value_(std::move(val))
        {}


    // Member Functions

        //- Return the value
        const Type& value() const noexcept
        {
            return value_;
        }

        //- Non-const access to the value
        Type& value() noexcept
        {
            return value_;
        }


    // Member Operators

        //- Implicit cast to the value
        operator const Type&() const noexcept
        {
            return value_;
        }

        //- Return the value
        const Type& operator[](const label) const noexcept
        {
            return value_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
