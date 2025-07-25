/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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
    attachMesh

Group
    grpMeshManipulationUtilities

Description
    Attach topologically detached mesh using prescribed mesh modifiers.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "attachPolyTopoChanger.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Attach topologically detached mesh using prescribed mesh modifiers"
    );

    #include "addOverwriteOption.H"

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const bool overwrite = args.found("overwrite");

    if (!overwrite)
    {
        ++runTime;
    }

    Info<< "Time = " << runTime.timeName() << nl
        << "Attaching sliding interface" << endl;

    attachPolyTopoChanger(mesh).attach();

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
