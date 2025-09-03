/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    kivaToFoam

Group
    grpMeshConversionUtilities

Description
    Convert a KIVA3v grid to OpenFOAM.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "Fstream.H"
#include "cellShape.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "oldCyclicPolyPatch.H"
#include "unitConversion.H"

using namespace Foam;

//- Supported KIVA versions
enum kivaVersions
{
    kiva3,
    kiva3v
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert a KIVA3v grid to OpenFOAM"
    );
    argList::noParallel();
    argList::addOption
    (
        "file",
        "name",
        "Specify alternative input file name - default is otape17"
    );
    argList::addOption
    (
        "version",
        "version",
        "Specify kiva version [kiva3|kiva3v] - default is '3v'"
    );
    argList::addOption
    (
        "zHeadMin",
        "scalar",
        "Minimum z-height for transferring liner faces to cylinder-head"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const fileName kivaFileName =
        args.getOrDefault<fileName>("file", "otape17");

    kivaVersions kivaVersion = kiva3v;
    if (args.found("version"))
    {
        const word versionName = args["version"];

        if (versionName == "kiva3")
        {
            kivaVersion = kiva3;
        }
        else if (versionName == "kiva3v")
        {
            kivaVersion = kiva3v;
        }
        else
        {
            FatalErrorInFunction
                << "KIVA file version " << versionName << " not supported"
                << exit(FatalError);

            args.printUsage();
            FatalError.exit(1);
        }
    }

    const scalar zHeadMin = args.getOrDefault<scalar>("zHeadMin", -GREAT);

    #include "readKivaGrid.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
