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
    Test-faMesh-try

Description
    Test for loading of faMesh

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "faMesh.H"
#include "polyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "addFaRegionOption.H"
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    #include "getFaRegionOption.H"

    autoPtr<faMesh> aMeshPtr = faMesh::TryNew(areaRegionName, mesh);

    Info<< "area-mesh: " << Switch::name(aMeshPtr) << nl;

    Info<< "\nEnd\n" << nl;

    return 0;
}


// ************************************************************************* //
