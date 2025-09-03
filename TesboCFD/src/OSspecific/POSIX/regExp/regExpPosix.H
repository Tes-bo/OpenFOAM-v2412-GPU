/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    Foam::regExpPosix

Description
    Wrapper around POSIX extended regular expressions
    with some additional prefix-handling. The prefix-handling is
    loosely oriented on PCRE regular expressions and provides a
    simple means of tuning the expressions.

    The prefixes are detected as \c (?...) at the beginning of
    the regular expression. Any unknown/unsupported prefixes are silently
    ignored.

    - "(?!i)" :
       one or more embedded pattern-match modifiers for the entire pattern.
    - the \c 'i' indicates ignore-case
    - the \c '!' (exclamation) indicates negated (inverted) matching
    .

SeeAlso
    The manpage regex(7) for more information about POSIX regular expressions.
    These differ somewhat from \c Perl and \c sed regular expressions.

SeeAlso
    Foam::regExp and Foam::regExpCxx

Warning
    This class should not be used directly.
    Use the Foam::regExp typedef instead.

\deprecated
    This class will be superseded by Foam::regExpCxx as compiler support
    for regular expressions continues to improve.

SourceFiles
    regExpPosix.C
    regExpPosixI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_regExpPosix_H
#define Foam_regExpPosix_H

#include "regExpCxx.H"
#include <regex.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class StringType> class SubStrings;

/*---------------------------------------------------------------------------*\
                         Class regExpPosix Declaration
\*---------------------------------------------------------------------------*/

class regExpPosix
{
    // Data Types

        //- Simple control types
        enum ctrlType { EMPTY = 0, NORMAL = 1, NEGATED = 2 };


    // Private Data

        //- Compiled regular expression
        regex_t* preg_;

        //- Track if input pattern is non-empty, negated etc.
        unsigned char ctrl_;


    // Private Member Functions

        //- Assign pattern
        bool set_pattern(const char* pattern, size_t len, bool ignoreCase);


public:

    // Public Types

        //- Type for matches - similar to std::smatch
        typedef SubStrings<std::string> results_type;


    // Static Member Data

        //- Grammar (unused) - for compatibility with Foam::regExpCxx
        static int grammar;


    // Static Member Functions

        //- Test if character is a regex meta-character
        inline static bool is_meta(const char c) noexcept
        {
            return regExpCxx::is_meta(c);
        }

        //- Test if string contains any (unquoted) meta-characters
        inline static bool is_meta
        (
            const std::string& str,
            const char quote = '\\'
        )
        {
            return regExpCxx::is_meta(str, quote);
        }


    // Public Classes

        //- Functor wrapper for testing meta-characters
        using meta = regExpCxx::meta;


    // Constructors

        //- Default construct
        inline regExpPosix() noexcept;

        //- Copy construct - disallowed
        regExpPosix(const regExpPosix&) = delete;

        //- Move construct
        inline regExpPosix(regExpPosix&& rgx) noexcept;

        //- Construct from character array, optionally ignore case
        inline explicit regExpPosix
        (
            const char* pattern,
            const bool ignoreCase = false
        );

        //- Construct from string, optionally ignore case
        inline explicit regExpPosix
        (
            const std::string& pattern,
            const bool ignoreCase = false
        );


    //- Destructor
    inline ~regExpPosix();


    // Member Functions

    // Access

        //- Return true if a precompiled expression does not exist
        inline bool empty() const noexcept;

        //- Return true if a precompiled expression exists
        inline bool exists() const noexcept;

        //- True if pattern matching is negated
        inline bool negated() const noexcept;

        //- Change pattern negation, return previous value
        inline bool negate(bool on) noexcept;

        //- The number of capture groups for a non-empty,
        //- non-negated expressions
        inline unsigned ngroups() const;


    // Editing

        //- Clear expression.
        //  \return True if expression had existed prior to the clear.
        bool clear();

        //- Swap contents
        inline void swap(regExpPosix& rgx);

        //- Compile pattern into a regular expression, optionally ignore case.
        //  \return True if the pattern was compiled
        inline bool set(const char* pattern, bool ignoreCase=false);

        //- Compile pattern into a regular expression, optionally ignore case.
        //  \return True if the pattern was compiled
        inline bool set(const std::string& pattern, bool ignoreCase=false);


    // Matching/Searching

        //- Find position within the text.
        //  \return The index where it begins or string::npos if not found
        //
        //  \note does not properly work with negated regex!
        std::string::size_type find(const std::string& text) const;

        //- True if the regex matches the entire text.
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        bool match(const std::string& text) const;

        //- True if the regex matches the text, set the matches.
        //  The first group starts at index 1 (0 is the entire match).
        //  The begin-of-line (^) and end-of-line ($) anchors are implicit
        //
        //  \note does not properly work with negated regex!
        bool match(const std::string& text, results_type& matches) const;

        //- Return true if the regex was found within the text
        inline bool search(const std::string& text) const;


    // Member Operators

        //- Perform match on text
        inline bool operator()(const std::string& text) const;

        //- Copy assignment - disallowed
        void operator=(const regExpPosix&) = delete;

        //- Move assignment
        inline void operator=(regExpPosix&& rgx);

        //- Assign and compile pattern from a character array.
        //  Matching is case sensitive.
        inline void operator=(const char* pattern);

        //- Assign and compile pattern from string.
        //  Matching is case sensitive.
        inline void operator=(const std::string& pattern);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "regExpPosixI.H"

#endif

// ************************************************************************* //
