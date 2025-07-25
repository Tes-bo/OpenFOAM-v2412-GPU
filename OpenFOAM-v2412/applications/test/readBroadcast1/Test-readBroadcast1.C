/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Test file reading with broadcast

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Fstream.H"
#include "Pstream.H"
#include "SpanStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noFunctionObjects();
    argList::noCheckProcessorDirectories();
    argList::addVerboseOption("additional information");
    argList::addBoolOption("fail", "fail if file cannot be opened");
    argList::addBoolOption("no-broadcast", "suppress broadcast contents");

    argList::addNote("Test master-only reading (with broadcast)");

    argList::addArgument("srcFile");

    #include "setRootCase.H"

    const bool syncPar = (UPstream::parRun() && !args.found("no-broadcast"));
    const bool optFail = args.found("fail");

    auto srcName = args.get<fileName>(1);

    if (srcName.has_ext("gz"))
    {
        srcName.remove_ext();
        Info<< "stripping extraneous .gz ending" << endl;
    }

    ICharStream is;

    {
        DynamicList<char> buffer;
        if (UPstream::master() || !syncPar)
        {
            if (optFail)
            {
                IFstream ifs(srcName, IOstreamOption::BINARY);

                if (!ifs.good())
                {
                    FatalIOErrorInFunction(srcName)
                        << "Cannot open file " << srcName
                        << exit(FatalIOError);
                }

                buffer = IFstream::readContents(ifs);
            }
            else
            {
                buffer = IFstream::readContents(srcName);
            }
        }

        if (syncPar)
        {
            // Prefer two broadcasts instead of serialize/de-serialize
            Pstream::broadcastList(buffer);
        }

        is.swap(buffer);
    }

    Pout<< "input:" << is.capacity() << endl;

    for (string line; is.getLine(line); /*nil*/)
    {
        Pout<< "L:" << is.lineNumber() << ": " << line.c_str() << nl;
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
