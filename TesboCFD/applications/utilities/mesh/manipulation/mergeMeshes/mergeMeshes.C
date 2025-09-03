/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    mergeMeshes

Group
    grpMeshManipulationUtilities

Description
    Merges two meshes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "mergePolyMesh.H"
#include "topoSet.H"
#include "processorMeshes.H"

using namespace Foam;

void getRootCase(fileName& casePath)
{
    casePath.clean();  // Remove unneeded ".."

    if (casePath.empty() || casePath == ".")
    {
        // handle degenerate form and '.'
        casePath = cwd();
    }
    else if (casePath[0] != '/' && casePath.name() == "..")
    {
        // avoid relative cases ending in '..' - makes for very ugly names
        casePath = cwd()/casePath;
        casePath.clean();  // Remove unneeded ".."
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge two meshes"
    );

    #include "addOverwriteOption.H"

    argList::addArgument("masterCase");
    argList::addOption
    (
        "masterRegion",
        "name",
        "Specify alternative mesh region for the master mesh"
    );

    argList::addArgument("addCase");
    argList::addOption
    (
        "addRegion",
        "name",
        "Specify alternative mesh region for the additional mesh"
    );
    argList::addOption
    (
        "resultTime",
        "time",
        "Specify a time for the resulting mesh"
    );

    argList args(argc, argv);
    if (!args.check())
    {
         FatalError.exit();
    }

    const bool overwrite = args.found("overwrite");

    auto masterCase = args.get<fileName>(1);
    auto addCase = args.get<fileName>(2);

    const word masterRegion =
        args.getOrDefault<word>("masterRegion", polyMesh::defaultRegion);

    const word addRegion =
        args.getOrDefault<word>("addRegion", polyMesh::defaultRegion);

    // Since we don't use argList processor directory detection, add it to
    // the casename ourselves so it triggers the logic inside TimePath.
    const fileName& cName = args.caseName();
    const auto pos = cName.find("processor");
    if (pos != string::npos && pos != 0)
    {
        fileName processorName = cName.substr(pos);
        masterCase += '/' + processorName;
        addCase += '/' + processorName;
    }


    getRootCase(masterCase);
    getRootCase(addCase);

    Info<< "Master:      " << masterCase << "  region " << masterRegion << nl
        << "mesh to add: " << addCase    << "  region " << addRegion << endl;

    #include "createTimes.H"

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName() << nl;

    Info<< "Create mesh\n" << endl;
    mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );

    Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName() << nl;
    Info<< "Create mesh\n" << endl;
    polyMesh meshToAdd
    (
        IOobject
        (
            addRegion,
            runTimeToAdd.timeName(),
            runTimeToAdd
        )
    );

    word meshInstance = masterMesh.pointsInstance();

    const bool specifiedInstance =
    (
        !overwrite
     && args.readIfPresent("resultTime", meshInstance)
    );

    if (specifiedInstance)
    {
        runTimeMaster.setTime(instant(meshInstance), 0);
    }
    else if (!overwrite)
    {
        runTimeMaster++;
    }

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();

    if (overwrite || specifiedInstance)
    {
        masterMesh.setInstance(meshInstance);
    }

    masterMesh.write();
    topoSet::removeFiles(masterMesh);
    processorMeshes::removeFiles(masterMesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
