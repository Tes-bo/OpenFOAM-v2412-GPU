/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
    Test-etcFiles

Description
    Test etcFiles functionality.
    Similar to foamEtcFile script, but automatically prunes nonexistent
    directories from the list.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "etcFiles.H"
#include "foamVersion.H"

using namespace Foam;

void printList(const fileNameList& list)
{
    for (const fileName& f : list)
    {
        Info<< f.c_str() << nl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noBanner();
    argList::noParallel();
    argList::noFunctionObjects();
    argList::removeOption("case");

    argList::addBoolOption
    (
        "config",
        "Print compile-time configuration values"
    );
    argList::addBoolOption
    (
        "all",
        "Return all files (otherwise stop after the first match)"
    );
    argList::addBoolOption
    (
        "list",
        "List directories or files to be checked"
    );
    argList::addBoolOption
    (
        "list-all",
        "List all directories (including non-existent ones)"
    );
    argList::addArgument("file...");

    argList::addNote
    (
        "Locate user/group/other file with semantics similar to the "
        "<etc>/fileName expansion."
    );

    argList args(argc, argv, false, true);

    // First handle no parameters
    if (args.size() == 1)
    {
        if (args.found("config"))
        {
            Info<<"config:project=" << foamVersion::configuredProjectDir << nl;
            Info<<"config:etc=" << foamVersion::configuredEtcDir << nl;
            return 0;
        }
        else if (args.found("list-all"))
        {
            fileNameList results = etcDirs(false);
            printList(results);
            return 0;
        }
        else if (args.found("list"))
        {
            fileNameList results = etcDirs();
            printList(results);
            return 0;
        }
        else
        {
            Info<<"Error: Missing filename" << endl;
            args.printUsage();
            return 1;
        }
    }


    // This should and will fail:
    //// Info<<"find bad file:" << nl
    ////     << findEtcFile("##BadName##", true) << "FAIL" << endl;


    const bool listAll = (args.found("all") || args.found("list"));

    int error = 0;

    for (int argi = 1; argi < args.size(); ++argi)
    {
        const std::string file = args[argi];
        fileNameList results = findEtcFiles(file);

        if (results.empty())
        {
            Info<<"Not found: " << file << nl;
            error = 2;
        }
        else if (listAll)
        {
            printList(results);
        }
        else
        {
            Info<<results[0].c_str() << nl;
        }
    }

    return error;
}


// ************************************************************************* //
