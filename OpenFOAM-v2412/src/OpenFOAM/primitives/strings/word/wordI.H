/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include <iostream>  // For std::cerr

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class PrimitiveType>
inline Foam::word Foam::word::printf
(
    const char* fmt,
    const PrimitiveType& val
)
{
    word output;
    string_printf(output, fmt, val);
    return output;
}


template<class PrimitiveType>
inline Foam::word Foam::word::printf
(
    const std::string& fmt,
    const PrimitiveType& val
)
{
    word output;
    string_printf(output, fmt, val);
    return output;
}


inline bool Foam::word::valid(char c)
{
    return
    (
        !isspace(c)
     && c != '"'   // string quote
     && c != '\''  // string quote
     && c != '/'   // path separator
     && c != ';'   // end statement
     && c != '{'   // beg block (eg, subdict)
     && c != '}'   // end block (eg, subdict)
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::word::word(const string& s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::word::word(string&& s, bool doStrip)
:
    string(std::move(s))
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::word::word(std::string&& s, bool doStrip)
:
    string(std::move(s))
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::word::word(const std::string& s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::word::word(const char* s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::word::word(const char* s, size_type len, bool doStrip)
:
    string(s, len)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline void Foam::word::stripInvalid()
{
    // Only strip when debug is active (potentially costly operation)
    if (debug && string::stripInvalid<word>(*this))
    {
        std::cerr
            << "word::stripInvalid() called for word "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::exit(1);
        }
    }
}


// FUTURE?
// inline Foam::word Foam::word::stem() const
// {
//     const auto i = find_ext();  // Or just rfind('.')
//
//     if (i == std::string::npos)
//     {
//         return *this;
//     }
//
//     return substr(0, i);
// }


inline Foam::word Foam::word::ext() const
{
    return string::ext();
}


inline Foam::word& Foam::word::ext(const word& ending)
{
    string::ext(ending);
    return *this;
}


inline Foam::word& Foam::word::replace_ext(const word& ending)
{
    string::remove_ext();
    string::ext(ending);
    return *this;
}


inline Foam::word Foam::word::lessExt() const
{
    const auto i = find_ext();  // Or just rfind('.')

    if (i == std::string::npos)
    {
        return *this;
    }

    return substr(0, i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::word& Foam::word::operator=(const word& s)
{
    // Self-assignment is a no-op
    if (this != &s)
    {
        assign(s);
    }
    return *this;
}


inline Foam::word& Foam::word::operator=(word&& s)
{
    // Self-assignment is a no-op
    if (this != &s)
    {
        assign(std::move(s));
    }
    return *this;
}


inline Foam::word& Foam::word::operator=(const string& s)
{
    assign(s);
    stripInvalid();
    return *this;
}


inline Foam::word& Foam::word::operator=(string&& s)
{
    assign(std::move(s));
    stripInvalid();
    return *this;
}


inline Foam::word& Foam::word::operator=(const std::string& s)
{
    assign(s);
    stripInvalid();
    return *this;
}


inline Foam::word& Foam::word::operator=(std::string&& s)
{
    assign(std::move(s));
    stripInvalid();
    return *this;
}


inline Foam::word& Foam::word::operator=(const char* s)
{
    assign(s);
    stripInvalid();
    return *this;
}


// ************************************************************************* //
