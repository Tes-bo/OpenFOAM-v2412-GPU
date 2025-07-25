/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    moveMesh

Group
    grpMeshManipulationUtilities

Description
    A solver utility for moving meshes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "motionSolver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "A solver utility for moving meshes"
    );

    argList::addOption
    (
        "deltaT",
        "time",
        "Override deltaT (eg, for accelerated motion)"
    );

    argList::addOption
    (
        "endTime",
        "time",
        "Override endTime (eg, for shorter tests)"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    scalar timeVal = 0;
    if (args.readIfPresent("deltaT", timeVal))
    {
        runTime.setDeltaT(timeVal);
    }

    if (args.readIfPresent("endTime", timeVal))
    {
        runTime.stopAt(Time::stopAtControls::saEndTime);
        runTime.setEndTime(timeVal);
    }

    autoPtr<motionSolver> motionPtr = motionSolver::New(mesh);

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        mesh.movePoints(motionPtr->newPoints());

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
