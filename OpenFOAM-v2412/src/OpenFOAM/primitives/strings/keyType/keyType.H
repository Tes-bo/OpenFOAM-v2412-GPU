/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::keyType

Description
    A class for handling keywords in dictionaries.

    A keyType is the keyword of a dictionary.
    It differs from word in that it also accepts patterns (regular expressions).
    It is very similar to wordRe, but doesn't store a regular expression.

SourceFiles
    keyType.C
    keyTypeI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_keyType_H
#define Foam_keyType_H

#include "word.H"
#include "wordRe.H"
#include "stdFoam.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class keyType;
class token;
Istream& operator>>(Istream& is, keyType& val);
Ostream& operator<<(Ostream& os, const keyType& val);

//- Hashing for keyType
template<> struct Hash<keyType> : string::hasher {};


/*---------------------------------------------------------------------------*\
                           Class keyType Declaration
\*---------------------------------------------------------------------------*/

class keyType
:
    public word
{
public:

    // Public Data Types

        //- Enumeration for the data type and search/match modes (bitmask)
        //  eg, (keyType::REGEX | keyType::RECURSIVE)
        enum option : unsigned char
        {
            // Base stored types
            LITERAL = 0,        //!< String literal
            REGEX   = 1,        //!< Regular expression

            // Variants for search/match only
            RECURSIVE = 0x80,   //!< Recursive search (eg, in dictionary)
            LITERAL_RECURSIVE = (LITERAL | RECURSIVE),
            REGEX_RECURSIVE   = (REGEX | RECURSIVE)
        };


private:

    // Private Data

        //- Treat keyType as literal, regex etc.
        //  Never contains RECURSIVE values.
        option type_;


public:

    // Static Data Members

        //- An empty keyType
        static const keyType null;


    // Constructors

        //- Default construct, empty literal
        inline keyType();

        //- Copy construct
        keyType(const keyType&) = default;

        //- Move construct
        keyType(keyType&&) = default;

        //- Implicit copy construct from word, treat as LITERAL
        inline keyType(const word& str);

        //- Implicit move construct from word, treat as LITERAL
        inline keyType(word&& str);

        //- Implicit copy construct from Foam::string, treat as REGEX
        inline keyType(const string& str);

        //- Implicit move construct from Foam::string, treat as REGEX
        inline keyType(string&& str);

        //- Copy construct from std::string with specified treatment
        inline keyType(const std::string& str, option opt);

        //- Move construct from std::string with specified treatment
        inline keyType(std::string&& str, option opt);

        //- Implicit construct from character array,
        //- with specified compile option (default is LITERAL)
        inline keyType(const char* str, option opt = option::LITERAL);

        //- Construct from Istream by reading a token
        //  Treat as regular expression if surrounded by quotation marks.
        explicit keyType(Istream& is);


    // Member Functions

        //- Test for valid keyType character?
        //  Like Foam::word, but with brace-brackets,
        //  which are valid for some regexs.
        inline static bool valid(const char c);


    // Access

        //- The keyType is treated as literal, not as pattern.
        inline bool isLiteral() const noexcept;

        //- The keyType is treated as a pattern, not as literal string.
        inline bool isPattern() const noexcept;


    // Infrastructure

        //- Inherit all regular string assign() methods
        using word::assign;

        //- Assign from word or string token.
        //  Words are treated as literals, strings as regex
        //  \return false if the token was the incorrect type
        bool assign(const token& tok);

        //- Change the representation, optionally stripping invalid word
        //- characters when changing to a literal
        inline void setType(option opt, bool adjust = false);

        //- Mark as regular expression
        inline bool compile() noexcept;

        //- Mark as literal string
        inline void uncompile() noexcept;

        //- Mark as literal string, optionally strip invalid word
        //- characters when changing to a literal
        inline void uncompile(bool adjust);


    // Editing

        //- Clear string and set as literal
        inline void clear();

        //- Swap contents. Self-swapping is a no-op.
        void swap(keyType& rhs);


    // Matching/Searching

        //- Smart match as regular expression or as a string.
        //  Optionally force a literal match only
        bool match(const std::string& text, bool literal=false) const;


    // Member Operators

        //- Perform smart match on text, as per match()
        //  Allows use as a predicate.
        inline bool operator()(const std::string& text) const;

        //- No assignment where type could be indeterminate
        void operator=(const std::string&) = delete;

        //- Copy assignment, retaining type (literal or regex)
        //  Self-assignment is a no-op.
        inline void operator=(const keyType& str);

        //- Move assignment, retaining type (literal or regex)
        //  Self-assignment is a no-op.
        inline void operator=(keyType&& str);

        //- Assign from word, treat as literal
        inline void operator=(const word& str);

        //- Assign from Foam::string, treat as regular expression
        inline void operator=(const string& str);

        //- Assign from character array, treat as literal
        inline void operator=(const char* str);


    // Housekeeping

        //- Deprecated(2019-08) construct as literal/regex
        //  \deprecated(2019-08) - use Construct with option
        FOAM_DEPRECATED_FOR(2019-08, "Construct with option")
        keyType(const std::string& s, bool isPattern)
        :
            keyType(s, (isPattern ? option::REGEX : option::LITERAL))
        {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, keyType& val);

//- Write operator
Ostream& operator<<(Ostream& os, const keyType& val);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "keyTypeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
