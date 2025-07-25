/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2015 OpenFOAM Foundation
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
    createExternalCoupledPatchGeometry.

Group
    grpPreProcessingUtilities

Description
    Generate the patch geometry (points and faces) for use
    with the externalCoupled functionObject.

Usage
    \verbatim
    createExternalCoupledPatchGeometry \<patchGroup\> [OPTION]
    \endverbatim

    \param -commsDir \<commsDir\> \n
    Specify an alternative communications directory (default is comms
    in the case directory)

    \param -region \<name\> \n
    Specify an alternative mesh region.

    \param -regions (\<name1\> .. \<nameN\>) \n
    Specify alternative mesh regions. The region names will be sorted
    alphabetically and a single composite name will be created
        \<nameX\>_\<nameY\>.._\<nameZ\>

    On execution, the combined patch geometry (points and faces) are output
    to the communications directory.

Note:
    The addressing is patch-local, i.e. point indices for each patch point
    used for face addressing starts at index 0.

See also
    functionObjects::externalCoupled

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "externalCoupled.H"
#include "regionProperties.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate the patch geometry (points and faces) for use"
        " with the externalCoupled functionObject."
    );

    argList::addOption
    (
        "region",
        "name",
        "Specify alternative mesh region"
    );
    argList::addOption
    (
        "regions",
        "(name1 .. nameN)",
        "Specify alternative mesh regions"
    );

    argList::addArgument("patchGroup");
    argList::addOption
    (
        "commsDir",
        "dir",
        "Specify communications directory (default is 'comms')"
    );
    #include "setRootCase.H"
    #include "createTime.H"

    wordList regionNames(1, polyMesh::defaultRegion);
    if (!args.readIfPresent("region", regionNames.first()))
    {
        args.readIfPresent("regions", regionNames);
    }

    const wordRe patchGroup(args.get<wordRe>(1));

    fileName commsDir(runTime.path()/"comms");
    args.readIfPresent("commsDir", commsDir);


    // Make sure region names are in canonical order
    stableSort(regionNames);


    PtrList<const fvMesh> meshes(regionNames.size());
    forAll(regionNames, regioni)
    {
        const word& regionName = regionNames[regioni];

        Info<< "Create mesh " << regionName
            << " for time = "
            << runTime.timeName() << nl << endl;

        meshes.set
        (
            regioni,
            new fvMesh
            (
                Foam::IOobject
                (
                    regionName,
                    runTime.timeName(),
                    runTime,
                    Foam::IOobject::MUST_READ
                )
            )
        );
    }


    functionObjects::externalCoupled::writeGeometry
    (
        UPtrList<const fvMesh>(meshes),
        commsDir,
        patchGroup
    );

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
