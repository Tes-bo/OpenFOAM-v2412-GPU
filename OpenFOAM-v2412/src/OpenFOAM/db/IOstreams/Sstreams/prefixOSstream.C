/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "prefixOSstream.H"
#include "Pstream.H"
#include "token.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

inline void Foam::prefixOSstream::checkWritePrefix()
{
    if (printPrefix_ && !prefix_.empty())
    {
        OSstream::write(prefix_.c_str());
        printPrefix_ = false;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::prefixOSstream::prefixOSstream
(
    std::ostream& os,
    const string& streamName,
    IOstreamOption streamOpt
)
:
    OSstream(os, streamName, streamOpt),
    printPrefix_(true),
    prefix_()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::prefixOSstream::print(Ostream& os) const
{
    os  << "prefixOSstream ";
    OSstream::print(os);
}


bool Foam::prefixOSstream::write(const token& tok)
{
    return OSstream::write(tok);
}


Foam::Ostream& Foam::prefixOSstream::write(const char c)
{
    checkWritePrefix();
    OSstream::write(c);

    // Reset prefix state on newline
    if (c == token::NL) printPrefix_ = true;
    return *this;
}


Foam::Ostream& Foam::prefixOSstream::writeQuoted
(
    const char* str,
    std::streamsize len,
    const bool quoted
)
{
    checkWritePrefix();
    OSstream::writeQuoted(str, len, quoted);
    return *this;
}


Foam::Ostream& Foam::prefixOSstream::write(const char* str)
{
    checkWritePrefix();
    OSstream::write(str);

    const size_t len = strlen(str);
    if (len > 0 && str[len-1] == token::NL)
    {
        // Reset prefix state on newline
        printPrefix_ = true;
    }

    return *this;
}


Foam::Ostream& Foam::prefixOSstream::write(const word& val)
{
    // Unquoted, and no newlines expected.
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const std::string& str)
{
    return writeQuoted(str.data(), str.size(), true);
}


Foam::Ostream& Foam::prefixOSstream::write(const int32_t val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const int64_t val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const float val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write(const double val)
{
    checkWritePrefix();
    return OSstream::write(val);
}


Foam::Ostream& Foam::prefixOSstream::write
(
    const char* buf,
    std::streamsize count
)
{
    checkWritePrefix();
    return OSstream::write(buf, count);
}


void Foam::prefixOSstream::indent()
{
    checkWritePrefix();
    OSstream::indent();
}


// ************************************************************************* //
