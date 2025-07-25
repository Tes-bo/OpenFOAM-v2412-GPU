/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

Application
    Test-SHA1

Description

\*---------------------------------------------------------------------------*/

#include "OSHA1stream.H"
#include "SpanStream.H"
#include "StringStream.H"
#include "dictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char * argv[])
{
    SHA1 sha;
    OSHA1stream osha;
    SHA1Digest shaDig;

    const std::string str("The quick brown fox jumps over the lazy dog");
    Info<< shaDig << " : empty" << nl;
    Info<< SHA1(str) << " : " << str << nl;

    Info<< osha.digest() << " : empty" << nl;

    osha<< "";
    Info<< osha.digest() << " : still empty" << nl;

    osha<< std::string();
    Info<< osha.digest() << " : still empty" << nl;

    sha.append(str);
    Info<< sha << " : appended to empty" << nl;

    osha<< str;
    Info<< osha.digest() << " : output << to empty" << nl;

    sha.clear();
    sha.append(str);
    shaDig = sha;

    sha.append("\n");
    Info<< sha << " : with trailing newline" << nl;
    Info<< shaDig << nl;

    if (sha == shaDig)
    {
        Info<<"SHA1 digests are identical (unexpected)\n";
    }
    else
    {
        Info<<"SHA1 digests are different (expected)\n";
    }
    Info<<"lhs:" << sha << " rhs:" << shaDig << endl;

    // start over:
    sha.clear();
    sha.append(str);

    SHA1 sha_A = sha;

    sha.append("\n");

    Info<< "digest1: " << sha_A << nl;
    Info<< "digest2: " << sha << nl;

    // start over:
    sha.clear();
    sha.append("\"");
    sha.append(str);
    sha.append("\"");

    Info<< "digest3: " << sha << nl;

    // try the output buffer interface
    {
        OSHA1stream os;

        os  << str;
        Info<< os.digest() << endl;

        os  << str;
        Info<< os.digest() << endl;

        os.reset();
        os  << "The quick brown fox jumps over the lazy dog";
        Info<< os.digest() << endl;
    }

    {
        dictionary dict
        (
            ICharStream
            (
                "parent { Default_Boundary_Region { type zeroGradient; } }"
                "inlet_1 { value inlet_1; }"
                "inlet_2 { value inlet_2; }"
                "inlet_3 { value inlet_3; }"
                "\"inlet_.*\" { value XXX; }"
            ) ()
        );

        Info<< "dict:" << endl;
        dict.write(Info, false);

        dictionary dict2(dict);

        OSHA1stream os;
        dict.write(os, false);
        Info<< os.digest() << endl;

        Info<< dict2.digest() << endl;
    }


    return 0;
}
