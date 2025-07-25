/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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
    Test-faMeshesRegistry

Description
    Basic tests for faMeshesRegistry

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "faMesh.H"
#include "faMeshesRegistry.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createPolyMesh.H"

    Info<< "mesh 0: " << mesh.sortedNames() << nl;

    faMeshesRegistry& reg =
        const_cast<faMeshesRegistry&>(faMeshesRegistry::New(mesh));

    // faMeshesRegistry faReg = faMeshesRegistry(mesh);

    faMesh mesh1(mesh, Foam::zero{});
    faMesh mesh2("mesh2", mesh, Foam::zero{});

    reg.write();

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
