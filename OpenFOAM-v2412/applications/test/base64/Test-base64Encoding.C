/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenCFD Ltd.
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

Application
    Test-base64Encoding

Description
    Test base64 encoding layer.

    Target values generated with "base64 encode ..." in Google.
    A simple independent source for comparison.

\*---------------------------------------------------------------------------*/

#include "base64Layer.H"
#include "List.H"
#include "Pair.H"

#include <sstream>

using namespace Foam;

bool test(const Pair<string>& unit)
{
    const string& input = unit.first();
    const string& expected = unit.second();

    std::ostringstream os;

    base64Layer b64(os);
    b64.write(input.data(), input.size());
    b64.close();

    const string encoded = os.str();

    Info<< input << nl;

    if (encoded == expected)
    {
        Info<< "  encoded: " << encoded << " (OK)" << nl
            << endl;
        return true;
    }
    else
    {
        Info<< "  encoded: " << encoded << " (ERROR)" << nl
            << "  expected: " << expected << nl
            << endl;

        return false;
    }
}


bool test(std::initializer_list<Pair<string>> list)
{
    bool good = true;

    for (const Pair<string>& t : list)
    {
        good = test(t) && good;
    }

    return good;
}


bool test(const UList<Pair<string>>& list)
{
    bool good = true;

    for (const Pair<string>& t : list)
    {
        good = test(t) && good;
    }

    return good;
}


void testMixed(std::ostream& os, const UList<Pair<string>>& list)
{
    base64Layer b64(os);

    os  << "<test-mixed>" << nl;

    int i=0;
    for (const Pair<string>& t : list)
    {
        const string& input = t.first();

        os  << "<input" << ++i << " value='" << input << "'>" << nl;
        os  << "  ";

        b64.write(input.data(), input.size());
        b64.close();

        os  << nl
            << "</input" << i << ">" << nl;
    }

    os  << "</test-mixed>" << nl
        <<  nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char * argv[])
{
    Info<< "Test base64 encode layer" << nl << endl;

    List<Pair<string>> testList
    {
        {
            "abcdef",  // 6 input, 8 output
            "YWJjZGVm"
        },
        {
            "OpenFOAM",
            "T3BlbkZPQU0="
        },
        {
            "OpenFOAM: The Open Source CFD Toolbox",
            "T3BlbkZPQU06IFRoZSBPcGVuIFNvdXJjZSBDRkQgVG9vbGJveA=="
        }
    };


    bool good = test(testList);


    // Test mixing output
    testMixed(std::cout, testList);

    if (good)
    {
        Info<< "All tests passed" << endl;
        return 0;
    }
    else
    {
        Info<< "One or more tests failed" << endl;
        return 1;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
