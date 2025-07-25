/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

Description

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "SpanStream.H"
#include "StringStream.H"

using namespace Foam;

template<class OStreamType>
void doTest()
{
    OStreamType os;
    os << "output with some values " << 1 << " entry" << endl;

    Info<< "contains:" << nl
        << os.str() << endl;
    os.rewind();

    Info<< "after rewind:" << nl
        << os.str() << endl;

    os << "####";

    Info<< "overwrite with short string:" << nl
        << os.str() << endl;

    os.reset();
    os << "%%%% reset";

    Info<< "after reset:" << nl
        << os.str() << endl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    Info<< "Begin test OStringStream" << endl;

    doTest<OStringStream>();
    // No reset() method:  doTest<OCharStream>();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
