/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SpanStream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

    #include "setRootCase.H"
    #include "createTime.H"

    Info<< "Create mesh, no clear-out\n" << endl;
    polyMesh mesh
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );


    Info<< "Cell centres" << nl << mesh.cellCentres() << endl;
    Info<< "Cell volumes" << nl << mesh.cellVolumes() << endl;
    Info<< "Cell shapes" << nl << mesh.cellShapes() << endl;
    Info<< "Cell face centres" << nl << mesh.faceCentres() << endl;

    // Test construct from cellShapes
    {
        pointField points(mesh.points());
        cellShapeList shapes(mesh.cellShapes());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        faceListList boundaryFaces(pbm.size());
        forAll(pbm, patchi)
        {
            boundaryFaces[patchi] = pbm[patchi];
        }
        wordList boundaryPatchNames(pbm.names());
        PtrList<dictionary> boundaryDicts(pbm.size());
        forAll(pbm, patchi)
        {
            OCharStream os;
            os << pbm[patchi];
            ISpanStream is(os.view());
            boundaryDicts.set(patchi, new dictionary(is));
        }

        word defaultBoundaryPatchName = "defaultFaces";
        word defaultBoundaryPatchType = emptyPolyPatch::typeName;

        polyMesh newMesh
        (
            IOobject
            (
                "newMesh",
                runTime.timeName(),
                runTime,
                Foam::IOobject::NO_READ
            ),
            std::move(points),
            shapes,
            boundaryFaces,
            boundaryPatchNames,
            boundaryDicts,
            defaultBoundaryPatchName,
            defaultBoundaryPatchType
        );

        Info<< "New cell centres" << nl << newMesh.cellCentres() << endl;
        Info<< "New cell volumes" << nl << newMesh.cellVolumes() << endl;
        Info<< "New cell shapes" << nl << newMesh.cellShapes() << endl;
        Info<< "New cell face centres" << nl << newMesh.faceCentres() << endl;
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
