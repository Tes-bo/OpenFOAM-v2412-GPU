/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    char

Description
    A character and a pointer to a character string.

SourceFiles
    char.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_char_H
#define Foam_primitives_char_H

#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Read single character
char readChar(Istream& is);

//- Read single character
Istream& operator>>(Istream& is, char& c);

//- Write single character
Ostream& operator<<(Ostream& os, const char c);

//- Write a nul-terminated C-string
Ostream& operator<<(Ostream& os, const char* str);

//- Test for whitespace (C-locale)
inline bool isspace(char c) noexcept
{
    return
    (
        c == ' '    // (0x20)  space (SPC)
     || c == '\t'   // (0x09)  horizontal tab (TAB)
     || c == '\n'   // (0x0a)  newline (LF)
     || c == '\v'   // (0x0b)  vertical tab (VT)
     || c == '\f'   // (0x0c)  feed (FF)
     || c == '\r'   // (0x0d)  carriage return (CR)
    );
}


/*---------------------------------------------------------------------------*\
                         Specialization pTraits<char>
\*---------------------------------------------------------------------------*/

//- Template specialisation for pTraits\<char\>
template<>
class pTraits<char>
{
    char p_;

public:

    // Static Data Members

        static const char* const typeName;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(char val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator char() const noexcept { return p_; }

        //- Access the value
        operator char&() noexcept { return p_; }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
