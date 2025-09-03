/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::ISstream::ISstream
(
    std::istream& is,
    const string& streamName,
    IOstreamOption streamOpt
)
:
    Istream(streamOpt),
    name_(streamName),
    is_(is)
{
    if (is_.good())
    {
        setOpened();
        setGood();
    }
    else
    {
        syncState();
    }
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

inline Foam::ISstream& Foam::ISstream::get(char& c)
{
    is_.get(c);
    syncState();

    if (c == '\n' && good())
    {
        ++lineNumber_;
    }

    return *this;
}


inline int Foam::ISstream::peek()
{
    return is_.peek();
}


inline Foam::ISstream& Foam::ISstream::getLine(std::string& str, char delim)
{
    std::getline(is_, str, delim);
    syncState();

    // Unlike with ignore(), cannot use gcount() to test success
    if (delim == '\n')
    {
        ++lineNumber_;
    }

    return *this;
}


inline std::streamsize Foam::ISstream::getLine(std::nullptr_t, char delim)
{
    is_.ignore(std::numeric_limits<std::streamsize>::max(), delim);
    std::streamsize count = is_.gcount();
    syncState();

    if (delim == '\n' && count > 0)
    {
        ++lineNumber_;
    }

    return count;
}


inline Foam::ISstream& Foam::ISstream::putback(const char c)
{
    if (c == '\n')
    {
        --lineNumber_;
    }

    if (!is_.putback(c))
    {
        setBad();
    }

    syncState();

    return *this;
}


// ************************************************************************* //
