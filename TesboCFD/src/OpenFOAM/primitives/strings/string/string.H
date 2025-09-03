/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
    Foam::string

Description
    A class for handling character strings derived from std::string.

    Strings may contain any characters and therefore are delimited by quotes
    for IO : "any list of characters".

    Used as a base class for word and fileName.

See also
    Foam::findEtcFile() for information about the site/user OpenFOAM
    configuration directory

SourceFiles
    stringI.H
    string.C
    stringIO.C
    stringTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_string_H
#define Foam_string_H

#include "char.H"
#include "Hasher.H"
#include <cstdlib>
#include <cstring>
#include <string>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class string;
class word;
class wordRe;
class Istream;
class Ostream;

template<class T> struct Hash;

/*---------------------------------------------------------------------------*\
                           Class string Declaration
\*---------------------------------------------------------------------------*/

class string
:
    public std::string
{
protected:

    // Protected Member Functions

        //- Find position of a file extension dot, return npos on failure.
        //  A wrapped version of find_last_of("./") with additional logic.
        inline static std::string::size_type find_ext(const std::string& str);

        //- Find position of a file extension dot, return npos on failure.
        //  A wrapped version of find_last_of("./") with additional logic.
        inline std::string::size_type find_ext() const;

        //- A printf-style formatter for a primitive.
        template<class PrimitiveType>
        static std::string::size_type string_printf
        (
            std::string& output,
            const char* fmt,
            const PrimitiveType& val
        );

        //- A printf-style formatter for a primitive.
        template<class PrimitiveType>
        static std::string::size_type string_printf
        (
            std::string& output,
            const std::string& fmt,
            const PrimitiveType& val
        );

        //- Return file name extension (part after last .)
        word ext() const;

        //- Append a '.' and the ending.
        //  The '.' and ending will not be added when the ending is empty,
        //  or when the object was or ended with a '/'.
        //
        //  \return True if append occurred.
        bool ext(const word& ending);

        //- Return true if it has an extension or simply ends with a '.'
        inline bool has_ext() const;

        //- Return true if the extension is the same as the given ending.
        //  No proper nullptr protection.
        inline bool has_ext(const char* ending) const;

        //- Return true if the extension is the same as the given ending.
        inline bool has_ext(const std::string& ending) const;

        //- Return true if the extension matches the given ending.
        bool has_ext(const wordRe& ending) const;

        //- Remove leading path, return true if string changed.
        inline bool remove_path();

        //- Remove extension, return true if string changed.
        inline bool remove_ext();


public:

    // Public Classes

        //- Hashing functor for string and derived string classes
        struct hasher
        {
            unsigned operator()(const std::string& str, unsigned seed=0) const
            {
                return Foam::Hasher(str.data(), str.length(), seed);
            }
        };

        //- Deprecated hashing functor - use hasher
        //  \deprecated(2021-04) - use hasher
        struct hash : string::hasher {};


    // Static Data Members

        //- The type name "string"
        static const char* const typeName;

        //- The debug flag
        static int debug;

        //- An empty string
        static const string null;


    // Constructors

        //- Default construct
        string() = default;

        //- Copy construct from std::string
        inline string(const std::string& str);

        //- Move construct from std::string
        inline string(std::string&& str);

        //- Construct as copy of character array
        inline string(const char* str);

        //- Construct as copy with a maximum number of characters
        inline string(const char* str, const size_type len);

        //- Construct from a single character
        inline explicit string(const char c);

        //- Construct fill copies of a single character
        inline string(const size_type len, const char c);

        //- Construct from Istream
        explicit string(Istream& is);


    // Static Member Functions

        //- Avoid masking the normal std::string length() method
        using std::string::length;

        //- Length of the character sequence (with nullptr protection)
        static inline std::string::size_type length(const char* s)
        {
            return (s ? strlen(s) : 0);
        }

        //- The length of the string
        static inline std::string::size_type length(const std::string& s)
        {
            return s.size();
        }

        //- Does the string contain valid characters only?
        template<class StringType>
        static inline bool valid(const std::string& str);

        //- Strip invalid characters from the given string
        template<class StringType>
        static inline bool stripInvalid(std::string& str);

        //- Return a valid String from the given string
        template<class StringType>
        static inline StringType validate(const std::string& str);


    // Member Functions

        //- Test for equality.
        //  \return True when strings match literally.
        inline bool match(const std::string& text) const;

        //- Avoid masking the normal std::string replace
        using std::string::replace;

        //- Replace first occurrence of sub-string s1 with s2,
        //- beginning at pos
        string& replace
        (
            const std::string& s1,
            const std::string& s2,
            size_type pos = 0
        );

        //- Replace all occurrences of sub-string s1 with s2,
        //- beginning at pos in the string.
        //  A no-op if s1 is empty.
        string& replaceAll
        (
            const std::string& s1,
            const std::string& s2,
            size_type pos = 0
        );

        //- Replace any occurrence of s1 characters with c2,
        //- beginning at pos in the string.
        //  A no-op if s1 is empty.
        string& replaceAny
        (
            const std::string& s1,
            const char c2,
            size_type pos = 0
        );

        //- Inplace expand initial tags, tildes, and all occurrences of
        //- environment variables as per stringOps::expand
        //
        //  Any unknown entries are removed silently if allowEmpty is true
        //  \sa
        //  Foam::findEtcFile
        string& expand(const bool allowEmpty = false);

        //- Remove repeated characters
        //  \return True if string changed
        bool removeRepeated(const char character);

        //- Remove the given text from the start of the string.
        //  \return True if the removal occurred
        bool removeStart(const std::string& text);

        //- Remove leading character, unless string is a single character
        //  \return True if the removal occurred
        bool removeStart(const char c);

        //- Remove the given text from the end of the string.
        //  \return True if the removal occurred
        bool removeEnd(const std::string& text);

        //- Remove trailing character, unless string is a single character
        //  \return True if the removal occurred
        bool removeEnd(const char c);


    // Editing

        //- Swap contents. Self-swapping is a no-op.
        inline void swap(std::string& str);


    // Member Operators

        //- Test for equality. Allows use as a predicate.
        //  \return True when strings match literally.
        inline bool operator()(const std::string& text) const;


    // Housekeeping

        //- True if string contains given character (cf. C++23)
        bool contains(char c) const noexcept
        {
            return (find(c) != std::string::npos);
        }

        //- True if string contains given [string view] substring (cf. C++23)
        bool contains(const std::string& s) const noexcept
        {
            return (find(s) != std::string::npos);
        }

        //- True if string contains given substring (cf. C++23)
        bool contains(const char* s) const
        {
            return (find(s) != std::string::npos);
        }

        //- True if string starts with given character (cf. C++20)
        bool starts_with(char c) const
        {
            return (!empty() && front() == c);
        }

        //- True if string starts with given [string view] prefix (C++20)
        bool starts_with(const std::string& s) const
        {
            return (size() >= s.size() && !compare(0, s.size(), s));
        }

        //- True if string starts with given prefix (C++20)
        bool starts_with(const char* s) const
        {
            const auto len = strlen(s);
            return (size() >= len && !compare(0, len, s, len));
        }

        //- True if string ends with given character (cf. C++20)
        bool ends_with(char c) const
        {
            return (!empty() && back() == c);
        }

        //- True if string ends with given [string view] suffix (cf. C++20)
        bool ends_with(const std::string& s) const
        {
            return (size() >= s.size() && !compare(size()-s.size(), npos, s));
        }

        //- True if string ends with given suffix (cf. C++20)
        bool ends_with(const char* s) const
        {
            const auto len = strlen(s);
            return (size() >= len && !compare(size()-len, npos, s, len));
        }

        //- Count the number of occurrences of the specified character
        //- in the string
        //  Partially deprecated (NOV-2017) in favour of stringOps::count
        size_type count(const char c) const;

        //- Deprecated(2019-11)
        //  \deprecated(2019-11) use starts_with instead
        bool startsWith(const std::string& s) const { return starts_with(s); }

        //- Deprecated(2019-11)
        //  \deprecated(2019-11) use ends_with instead
        bool endsWith(const std::string& s) const { return ends_with(s); }

        //- Deprecated(2019-11)
        //  \deprecated(2019-11) use removeEnd instead
        bool removeTrailing(const char c) { return removeEnd(c); }
};


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

//- Hashing for Foam::string
template<> struct Hash<string> : string::hasher {};

//- Hashing for std:::string
template<> struct Hash<std::string> : string::hasher {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// IOstream Operators

//- Read operator
Istream& operator>>(Istream& is, string& val);

//- Write operator
Ostream& operator<<(Ostream& os, const string& val);

//- Write operator
Ostream& operator<<(Ostream& os, const std::string& val);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "stringI.H"

#ifdef NoRepository
    #include "stringTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
