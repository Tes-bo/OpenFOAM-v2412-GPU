/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

Primitive
    bool

Description
    System bool

SourceFiles
    bool.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_bool_H
#define Foam_primitives_bool_H

#include "pTraits.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// IOstream Operators

//- Read bool from stream using Foam::Switch(Istream&)
Istream& operator>>(Istream& is, bool& b);

//- Write bool to a stream as an label (integer) value
Ostream& operator<<(Ostream& os, const bool b);

//- Read bool from stream using Foam::Switch(Istream&)
bool readBool(Istream& is);


/*---------------------------------------------------------------------------*\
                         Specialization pTraits<bool>
\*---------------------------------------------------------------------------*/

//- Template specialisation for pTraits\<bool\>
template<>
class pTraits<bool>
{
    bool p_;

public:

    // Typedefs

        //- Component type
        typedef bool cmptType;

        //- Magnitude type
        typedef bool magType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of bool is 0
        static constexpr direction rank = 0;

        //- Number of components in bool is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const bool zero;
        static const bool one;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(bool val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator bool() const noexcept { return p_; }

        //- Access the value
        operator bool&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
