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
    foamCellZoneToVTK.C

Description
    Write tet-decomposed OpenFOAM mesh in VTK format.
    For diagnostic purposes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "foamVtkInternalMeshWriter.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write OpenFOAM cellZone mesh to VTK"
    );
    argList::addOption
    (
        "cellZone",
        "name",
        "Convert mesh subset corresponding to specified cellZone"
    );
    argList::addBoolOption
    (
        "list",
        "List names of cellZones and exit"
    );
    timeSelector::addOptions();

    #include "setRootCase.H"

    word cellZoneName;
    args.readIfPresent("cellZone", cellZoneName);

    const bool optList = args.found("list");

    if (optList)
    {
        if (!cellZoneName.empty())
        {
            Info<< "specify -list or -cellZone, but not both!" << nl;
            return 1;
        }
    }
    else if (cellZoneName.empty())
    {
        Info<< "Did not specify a cellZone!!" << nl;
        return 1;
    }

    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName("zonemesh-" + cellZoneName);

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

            if (optList)
            {
                Info<< "cellZones:" << nl;
                for (const word& name : mesh.cellZones().names())
                {
                    Info<< "    " << name << nl;
                }
            }
            else
            {
                const cellZone* zonePtr =
                    mesh.cellZones().cfindZone(cellZoneName);

                Info<< "cellZone " << cellZoneName;
                if (!zonePtr)
                {
                    Info<< " ... not found" << nl;
                    continue;
                }
                Info<< nl;

                const cellZone& zn = *zonePtr;

                // Define a subset
                vtk::vtuCells vtuCells;
                vtuCells.reset(mesh, zn);

                vtk::internalMeshWriter writer
                (
                    mesh,
                    vtuCells,
                    fileName
                    (
                        mesh.time().globalPath() / meshName
                    )
                );

                writer.writeGeometry();

                writer.beginCellData();
                writer.writeProcIDs();

                Info<< "Wrote " << writer.output().name() << nl;
            }
        }
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
