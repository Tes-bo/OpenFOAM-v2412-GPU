/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    foamMeshToTet-vtk

Description
    Write tet-decomposed OpenFOAM mesh in VTK format.
    For diagnostic purposes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"

namespace Foam
{
    void writeVTKtetMesh(const fileName& output, const polyMesh& mesh);
}

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write tet-decomposed OpenFOAM mesh in VTK"
    );
    argList::noParallel();
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = "tetmesh";
    if (args.found("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    #include "createPolyMesh.H"

    forAll(timeDirs, timei)
    {
        runTime.setTime(timeDirs[timei], timei);

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timei || state != polyMesh::UNCHANGED)
        {
            fileName meshName(exportName);
            if (state != polyMesh::UNCHANGED)
            {
                meshName += '_' + runTime.timeName();
            }

            writeVTKtetMesh(meshName, mesh);
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
