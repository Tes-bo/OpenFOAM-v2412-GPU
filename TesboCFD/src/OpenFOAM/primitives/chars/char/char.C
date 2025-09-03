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

\*---------------------------------------------------------------------------*/

#include "char.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<char>::typeName = "char";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pTraits<char>::pTraits(Istream& is)
{
    is >> p_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

char Foam::readChar(Istream& is)
{
    char c;
    is.read(c);
    return c;
}


Foam::Istream& Foam::operator>>(Istream& is, char& c)
{
    is.read(c);
    is.check(FUNCTION_NAME);
    return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char c)
{
    os.write(c);
    os.check(FUNCTION_NAME);
    return os;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char* str)
{
    os.write(str);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
