/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Test-fileOperation1

Description
   Test string parsing and other bits for fileOperation

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "fileOperation.H"
#include "SubList.H"
#include "IOobject.H"
#include "IOstreams.H"
#include "OSspecific.H"


using namespace Foam;

word toString(const fileOperation::procRangeType& group)
{
    if (group.empty())
    {
        return word::null;
    }
    return Foam::name(group.first()) + "-" + Foam::name(group.last());
}


void testSplitPath(const fileName& pathName)
{
    fileName path, procDir, local;
    fileOperation::procRangeType group;
    label nProcs;

    const label proci =
        fileOperation::splitProcessorPath
        (
            pathName,
            path,
            procDir,
            local,
            group,
            nProcs
        );


    Info<< nl
        << "Input    = " << pathName << nl
        << "  path   = " << path << nl
        << "  proc   = " << procDir << nl
        << "  local  = " << local << nl
        << "  group  = " << group << " = " << toString(group) << nl
        << "  proci   = " << proci << nl
        << "  nProcs  = " << nProcs << nl;
}


void testSplitPaths(std::initializer_list<const char* const> dirNames)
{
    for (const auto& dirName : dirNames)
    {
        testSplitPath(fileName(dirName));
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addArgument("fileName .. fileNameN");
    argList::addOption("istream", "file", "test Istream values");


    testSplitPaths
    ({
        "foo/bar",
        "foo/processor5/system",
        "foo/processors100_0-5/constant",
        "foo/processors20_12-16/constant",
        "/new-processor-gen/case1/processors20",
        "/new-processor-gen/case1/processors100_0-5/constant",
        "/new-processor-gen/case1/processors/input",
        "devel/processor/ideas/processor0/system",

        "/path/processor0Generation1/case1/processor10/input",

        "path/processors100_ab-cd/constant",
        "path/processors100_a11-d00/constant",
    });


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
