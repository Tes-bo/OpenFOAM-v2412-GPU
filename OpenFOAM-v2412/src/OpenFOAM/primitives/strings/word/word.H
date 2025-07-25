/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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
    Foam::word

Description
    A class for handling words, derived from Foam::string.

    A word is a string of characters without whitespace, quotes, slashes,
    semicolons or brace brackets. Words are delimited by whitespace.

SourceFiles
    wordI.H
    word.C
    wordIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_word_H
#define Foam_word_H

#include "string.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class word;
Istream& operator>>(Istream& is, word& val);
Ostream& operator<<(Ostream& os, const word& val);

//- Hashing for word
template<> struct Hash<word> : string::hasher {};


/*---------------------------------------------------------------------------*\
                           Class word Declaration
\*---------------------------------------------------------------------------*/

class word
:
    public string
{
public:

    // Static Data Members

        //- The typeName
        static const char* const typeName;

        //- Debugging
        static int debug;

        //- An empty word
        static const word null;


    // Constructors

        //- Default construct
        word() = default;

        //- Copy construct
        word(const word&) = default;

        //- Move construct
        word(word&& w) = default;

        //- Copy construct from Foam::string
        inline word(const string& s, bool doStrip=true);

        //- Move construct from Foam::string
        inline word(string&& s, bool doStrip=true);

        //- Copy construct from std::string
        inline word(const std::string& s, bool doStrip=true);

        //- Move construct from std::string
        inline word(std::string&& s, bool doStrip=true);

        //- Copy from character array
        inline word(const char* s, bool doStrip=true);

        //- Copy from buffer for a maximum number of characters
        inline word(const char* s, size_type len, bool doStrip);

        //- Construct from Istream
        explicit word(Istream& is);


    // Member Functions

        //- Use a printf-style formatter for a primitive.
        //  The representation is not checked for valid characters -
        //  it is assumed that the caller knows what they are doing
        template<class PrimitiveType>
        inline static word printf
        (
            const char* fmt,
            const PrimitiveType& val
        );

        //- Use a printf-style formatter for a primitive.
        //  The representation is not checked for valid characters -
        //  it is assumed that the caller knows what they are doing
        template<class PrimitiveType>
        inline static word printf
        (
            const std::string& fmt,
            const PrimitiveType& val
        );

        //- Is this character valid for a word?
        inline static bool valid(char c);

        //- Construct validated word (no invalid characters).
        //  Optionally prefix any leading digit with '_' to have words
        //  that work nicely as dictionary keywords.
        static word validate(const std::string& s, const bool prefix=false);

        //- Construct validated word (no invalid characters) from a sequence
        //- of characters in the range [first,last),
        //  Optionally prefix any leading digit with '_'.
        static word validate
        (
            const char* first,
            const char* last,
            const bool prefix=false
        );

        //- Strip invalid characters from this word
        //  Trips an abort on invalid characters for debug 2 or greater
        inline void stripInvalid();


    // File-like Functions

        //- Return file name extension (part after last .)
        inline word ext() const;

        //- Append a '.' and the ending, and return the object.
        //  The '.' and ending will not be added when the ending is empty,
        //  or when the file name is empty or ended with a '/'.
        inline word& ext(const word& ending);

        //- Remove extension (if any) and append a new one
        inline word& replace_ext(const word& ending);

        //- Return word without extension (part before last .)
        inline word lessExt() const;

        //- Remove extension, return true if string changed.
        using string::remove_ext;

        //- Various checks for extensions
        using string::has_ext;


    // Member Operators

    // Assignment

        //- Copy assignment, no character validation required.
        //  Self-assignment is a no-op
        inline word& operator=(const word& s);

        //- Move assignment, no character validation required
        //  Self-assignment is a no-op
        inline word& operator=(word&& s);

        //- Copy assignment from Foam::string, stripping invalid characters
        inline word& operator=(const string& s);

        //- Move assignment from Foam::string, stripping invalid characters
        inline word& operator=(string&& s);

        //- Copy assignment from std::string, stripping invalid characters
        inline word& operator=(const std::string& s);

        //- Move assignment from std::string, stripping invalid characters
        inline word& operator=(std::string&& s);

        //- Copy, stripping invalid characters
        inline word& operator=(const char* s);


    // Housekeeping

        //- Same as has_ext()
        bool hasExt() const { return has_ext(); }

        //- Same as has_ext()
        bool hasExt(const std::string& s) const { return has_ext(s); }

        //- Same as remove_ext()
        bool removeExt() { return remove_ext(); }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, word& val);

//- Write operator
Ostream& operator<<(Ostream& os, const word& val);


// * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * * //

//- Join words as camelCase, capitalizing the first letter of b.
//  No effect if either argument is empty.
word operator&(const word& a, const word& b);


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- A word representation of a memory address as hexadecimal.
//  No special handling of nullptr (renders as 0x0)
word name(const void* ptr);


//- Extract name (as a word) from an object, typically using its name() method.
template<class T>
struct nameOp
{
    word operator()(const T& obj) const
    {
        return obj.name();
    }

    //- Less-compare two objects by their name() method - for sorting
    bool operator()(const T& a, const T& b) const
    {
        return (a.name() < b.name());
    }
};


//- Extract type (as a word) from an object, typically using its type() method.
template<class T>
struct typeOp
{
    word operator()(const T& obj) const
    {
        return obj.type();
    }

    //- Less-compare two objects by their type() method - for sorting
    bool operator()(const T& a, const T& b) const
    {
        return (a.type() < b.type());
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
