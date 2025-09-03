/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "IOstream.H"
#include "error.H"
#include "argList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::fileName Foam::IOstream::staticName_("stream");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::IOstream::name() const
{
    return staticName_;
}


Foam::fileName Foam::IOstream::relativeName() const
{
    return argList::envRelativePath(this->name());
}


bool Foam::IOstream::check(const char* operation) const
{
    return fatalCheck(operation);
}


bool Foam::IOstream::fatalCheck(const char* operation) const
{
    const bool ok = !bad();

    if (!ok)
    {
        FatalIOErrorInFunction(*this)
            << "error in IOstream " << relativeName()
            << " for operation " << operation
            << exit(FatalIOError);
    }

    return ok;
}


void Foam::IOstream::print(Ostream& os) const
{
    os  << "IOstream: " << "Version "  << version() << ", format "
        << format() << ", line " << lineNumber();

    if (opened())
    {
        os  << ", OPENED";
    }

    if (closed())
    {
        os  << ", CLOSED";
    }

    if (good())
    {
        os  << ", GOOD";
    }

    if (eof())
    {
        os  << ", EOF";
    }

    if (fail())
    {
        os  << ", FAIL";
    }

    if (bad())
    {
        os  << ", BAD";
    }

    os  << endl;
}


void Foam::IOstream::print(Ostream& os, const int streamState) const
{
    if (streamState == std::ios_base::goodbit)
    {
        os  << "goodbit set : the last operation on stream succeeded" << endl;
    }
    else if (streamState & std::ios_base::badbit)
    {
        os  << "badbit set : characters possibly lost" << endl;
    }
    else if (streamState & std::ios_base::failbit)
    {
        os  << "failbit set : some type of formatting error" << endl;
    }
    else if (streamState & std::ios_base::eofbit)
    {
        os  << "eofbit set : at end of stream" << endl;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<IOstream>& iproxy
)
{
    (*iproxy).print(os);
    return os;
}


// ************************************************************************* //
