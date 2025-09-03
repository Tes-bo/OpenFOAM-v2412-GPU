/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::fileName::fileName(const word& s)
:
    string(s)
{}


inline Foam::fileName::fileName(word&& s)
:
    string(std::move(s))
{}


inline Foam::fileName::fileName(const string& s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(string&& s, bool doStrip)
:
    string(std::move(s))
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const std::string& s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(std::string&& s, bool doStrip)
:
    string(std::move(s))
{
    if (doStrip)
    {
        stripInvalid();
    }
}


inline Foam::fileName::fileName(const char* s, bool doStrip)
:
    string(s)
{
    if (doStrip)
    {
        stripInvalid();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::fileName::valid(char c)
{
    return
    (
        c != '"'    // string quote
     && c != '\''   // string quote
     && (!isspace(c) || (allowSpaceInFileName && c == ' '))
    );
}


inline void Foam::fileName::stripInvalid()
{
    // Only strip when debug is active (potentially costly operation)
    if (debug && string::stripInvalid<fileName>(*this))
    {
        std::cerr
            << "fileName::stripInvalid() called for invalid fileName "
            << this->c_str() << std::endl;

        if (debug > 1)
        {
            std::cerr
                << "    For debug level (= " << debug
                << ") > 1 this is considered fatal" << std::endl;
            std::exit(1);
        }

        removeRepeated('/');
        removeEnd('/');
    }
}


inline bool Foam::fileName::isAbsolute(const std::string& str)
{
    return !str.empty() &&
    (
        // Starts with '/', but also accept '\\' since it will be
        // converted to a generic '/' or it is part of a (windows)
        // UNC '\\server-name\path'
        // - accept even on non-windows systems

        (str[0] == '/' || str[0] == '\\')

#ifdef _WIN32
         // Filesytem root - eg, d:/path or d:\path
     || (
            (str.length() > 2 && str[1] == ':')
         && (str[2] == '/' || str[2] == '\\')
        )
#endif
    );
}


inline bool Foam::fileName::isAbsolute() const
{
    return fileName::isAbsolute(*this);
}


inline bool Foam::fileName::isBackup() const
{
    return isBackup(*this);
}


inline bool Foam::fileName::has_path() const
{
    return contains('/');
}


inline std::string Foam::fileName::path(const std::string& str)
{
    const auto i = str.rfind('/');

    if (i == std::string::npos)
    {
        return ".";
    }
    else if (i)
    {
        return str.substr(0, i);
    }

    return "/";
}


inline Foam::fileName Foam::fileName::path() const
{
    return fileName::path(*this);
}


inline std::string Foam::fileName::name(const std::string& str)
{
    const auto i = str.rfind('/');

    if (i == std::string::npos)
    {
        return str;
    }

    return str.substr(i+1);
}


inline Foam::word Foam::fileName::name() const
{
    return fileName::name(*this);
}


Foam::word Foam::fileName::ext() const
{
    return string::ext();
}


inline Foam::word Foam::fileName::stem() const
{
    return fileName::stem(*this);
}


inline Foam::fileName& Foam::fileName::ext(const word& ending)
{
    string::ext(ending);
    return *this;
}


inline Foam::fileName& Foam::fileName::replace_ext(const word& ending)
{
    string::remove_ext();
    string::ext(ending);
    return *this;
}


inline Foam::fileName Foam::fileName::lessExt() const
{
    const auto i = find_ext();

    if (i == std::string::npos)
    {
        return *this;
    }

    return substr(0, i);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline Foam::fileName& Foam::fileName::operator=(const fileName& str)
{
    // Self-assignment is a no-op
    if (this != &str)
    {
        assign(str);
    }
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(fileName&& str)
{
    // Self-assignment is a no-op
    if (this != &str)
    {
        assign(std::move(str));
    }
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(const word& str)
{
    assign(str);
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(word&& str)
{
    assign(std::move(str));
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(const string& str)
{
    assign(str);
    stripInvalid();
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(string&& str)
{
    assign(std::move(str));
    stripInvalid();
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(const std::string& str)
{
    assign(str);
    stripInvalid();
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(std::string&& str)
{
    assign(std::move(str));
    stripInvalid();
    return *this;
}


inline Foam::fileName& Foam::fileName::operator=(const char* str)
{
    assign(str);
    stripInvalid();
    return *this;
}


// ************************************************************************* //
