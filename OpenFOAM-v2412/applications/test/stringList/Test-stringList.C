/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "stringListOps.H"
#include "ListOps.H"
#include "FlatOutput.H"
#include "IOstreams.H"
#include "SpanStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    stringList strings
    {
        "hello",
        "heello",
        "heeello",
        "bye",
        "bbye",
        "bbbye",
        "okey",
        "okkey",
        "okkkey",
    };
    labelList matches;

    wordRes matcher1(ICharStream("( okey \"[hy]e+.*\" )")());

    Info<< "stringList " << strings << nl;

    {
        keyType key(".*ee.*", keyType::REGEX);
        matches = findMatchingStrings(regExp(key), strings);

        Info<< "matches found for regexp " << key << " :" << nl
            << matches << nl;

        for (const label idx : matches)
        {
            Info<< " -> " << strings[idx] << nl;
        }
    }

    Info<< "Match found using ListOps = "
        << ListOps::found(strings, regExp(".*ee.*")) << nl;

    Info<< "First index = "
        << ListOps::find(strings, regExp(".*ee.*")) << nl;

    Info<< endl;

    matches = findMatchingStrings(matcher1, strings);

    Info<< "matching " << flatOutput(matcher1) << " => "
        << matcher1.matching(strings) << nl;
    Info<< "matches found for " << flatOutput(matcher1) << " => "
        << matches << nl;

    for (const label idx : matches)
    {
        Info<< " -> " << strings[idx] << nl;
    }
    Info<< endl;

    stringList subLst = subsetStrings(regExp(".*ee.*"), strings);
    Info<< "subset stringList: " << subLst << nl;

    subLst = subsetStrings(matcher1, strings);
    Info<< "subset stringList: " << subLst << nl;

    inplaceSubsetStrings(matcher1, strings);
    Info<< "subsetted stringList: " << strings << nl;

    inplaceSubsetStrings(regExp(".*l.*"), strings);
    Info<< "subsetted stringList: " << strings << nl;

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
