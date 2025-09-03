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
    foamToFireMesh

Description
    Write an OpenFOAM mesh in AVL/FIRE fpma format

Usage
    \b foamToFireMesh [OPTION]

    Options:
      - \par -ascii
        Write in ASCII format instead of binary

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1 (ie, no scaling).

See also
    Foam::cellTable, Foam::meshWriter and Foam::fileFormats::FIREMeshWriter

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "FIREMeshWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write an OpenFOAM mesh in AVL/FIRE fpma format"
    );
    argList::noParallel();
    timeSelector::addOptions();

    argList::addBoolOption
    (
        "ascii",
        "Write in ASCII format instead of binary"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1 (none)"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = meshWriter::defaultMeshName;
    if (args.found("case"))
    {
        exportName += '-' + args.globalCaseName();
    }


    // write control options
    // ~~~~~~~~~~~~~~~~~~~~~
    fileFormats::FIREMeshWriter::binary = !args.found("ascii");

    // Default: no rescaling
    scalar scaleFactor = 1;
    if (args.readIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }

    #include "createPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            fileFormats::FIREMeshWriter writer(mesh, scaleFactor);

            fileName meshName(exportName);
            if (state != polyMesh::UNCHANGED)
            {
                meshName += '_' + runTime.timeName();
            }

            writer.write(meshName);
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
