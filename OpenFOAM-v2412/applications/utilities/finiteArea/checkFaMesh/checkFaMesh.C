/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    makeFaMesh

Description
    Check a finiteArea mesh

Original Authors
    Zeljko Tukovic, FAMENA
    Hrvoje Jasak, Wikki Ltd.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "argList.H"
#include "faMesh.H"
#include "polyMesh.H"
#include "areaFaMesh.H"
#include "edgeFaMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "processorFaPatch.H"
#include "foamVtkIndPatchWriter.H"
#include "foamVtkLineWriter.H"
#include "faMeshTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Check a finiteArea mesh"
    );

    argList::addBoolOption
    (
        "write-vtk",
        "Write mesh as a vtp (vtk) file for display or debugging"
    );

    argList::addOption
    (
        "geometryOrder",
        "N",
        "Test different geometry order - experimental!!",
        true  // Advanced option
    );

    #include "addRegionOption.H"
    #include "addFaRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    int geometryOrder(1);
    if (args.readIfPresent("geometryOrder", geometryOrder))
    {
        Info<< "Setting faMesh::geometryOrder = " << geometryOrder << nl
            << "(experimental)" << nl << endl;

        faMesh::geometryOrder(geometryOrder);
    }

    #include "createNamedFaMesh.H"

    Info<< "Time = " << runTime.timeName() << nl << endl;

    // Mesh information (verbose)
    faMeshTools::printMeshChecks(aMesh);

    if (args.found("write-vtk"))
    {
        #include "faMeshWriteVTK.H"
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
