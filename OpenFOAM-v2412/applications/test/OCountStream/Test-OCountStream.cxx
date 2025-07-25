/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Test null and counting output streams

\*---------------------------------------------------------------------------*/

#include "OCountStream.H"
#include "SpanStream.H"
#include "StringStream.H"
#include "Fstream.H"
#include "IOstreams.H"
#include "argList.H"

using namespace Foam;

template<class OS>
void generateOutput(OS& os)
{
    for (label i = 0; i < 50; ++i)
    {
        os  << 1002 << " " << "abcd" << " "
            << "def" << " " << 3.14159 << ";\n";
    }
}


void printInfo(OSstream& os)
{
    Info<< "name: " << os.name() << " : " << os.stdStream().tellp() << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addOption("write", "file", "test writing to file");

    #include "setRootCase.H"

    OCountStream cnt;
    OCharStream cstr;
    OStringStream sstr;
    ocountstream plain;

    generateOutput(cstr);
    generateOutput(sstr);
    generateOutput(cnt);
    generateOutput(plain);

    cnt.print(Info);

    Info<< "counter state: " << (cnt.stdStream().rdstate()) << nl
        << "via char-stream: " << label(cstr.view().size()) << " chars" << nl
        << "via string-stream: " << label(sstr.count()) << " chars" << nl
        << "via ocountstream: " << plain.count() << " chars" << endl;

    fileName outputName;
    args.readIfPresent("write", outputName);

    if (outputName.size())
    {
        IOstreamOption streamOpt;

        if (outputName.has_ext("gz"))
        {
            outputName.remove_ext();
            streamOpt.compression(IOstreamOption::COMPRESSED);
        }


        OFstream os1(outputName, streamOpt);
        OFstream os2(nullptr);   // A /dev/null equivalent
        OFstream os3("/dev/null");

        // Doubled output
        generateOutput(os1); generateOutput(os1);
        generateOutput(os2); generateOutput(os2);
        generateOutput(os3); generateOutput(os3);

        Info<< nl
            << "doubled output" << nl;

        printInfo(os1);
        printInfo(os2);
        printInfo(os3);

        // Rewind and test single output
        os1.rewind(); generateOutput(os1);
        os2.rewind(); generateOutput(os2);
        os3.rewind(); generateOutput(os3);

        Info<< nl
            << "single output" << nl;

        printInfo(os1);
        printInfo(os2);
        printInfo(os3);
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
