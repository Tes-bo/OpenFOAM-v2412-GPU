/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    slicedFieldTest

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "SlicedGeometricField.H"
#include "slicedFvPatchFields.H"
#include "slicedSurfaceFields.H"
#include "slicedVolFields.H"

#include "areaFaMesh.H"
#include "areaFields.H"
#include "edgeFields.H"
#include "slicedAreaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addBoolOption
    (
        "finite-area",
        "Test finite-area mesh/fields"
    );
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"

    if (args.found("finite-area"))
    {
        autoPtr<faMesh> faMeshPtr(faMesh::TryNew(mesh));

        if (!faMeshPtr)
        {
            Info<< "Stop: no finiteArea" << nl;
            return 0;
        }

        auto& aMesh = faMeshPtr();

        Info<< "Test with finiteArea" << nl;

        (void)aMesh.areaCentres(),
        (void)aMesh.faceAreaNormals();

        vectorField flatBoundary(aMesh.nBoundaryEdges(), Zero);

        {
            const auto& bfld = aMesh.faceAreaNormals().boundaryField();

            forAll(aMesh.boundary(), patchi)
            {
                const auto& src = bfld[patchi];
                const auto& p = aMesh.boundary()[patchi];

                vectorList::subList out = p.boundarySlice(flatBoundary);

                if (out.size() == src.size())
                {
                    out = src;
                }
            }
        }
        flatBoundary *= 100;


        slicedAreaVectorField foo
        (
            IOobject
            (
                "centres",
                runTime.timeName(),
                aMesh.thisDb(),
                IOobject::NO_REGISTER
            ),
            aMesh,
            dimLength,
            aMesh.areaCentres(),
            flatBoundary
        );

        Info<< "Weird combination of centres and normals!" << nl << nl;
        foo.writeData(Info.stream());

        Info<< "Done" << nl;
        return 0;
    }

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::MUST_READ
        ),
        mesh
    );

    slicedVolVectorField C
    (
        IOobject
        (
            "C2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimLength,
        mesh.cellCentres(),
        mesh.faceCentres()
    );

    Info<< C << endl;
    Info<< (C & U) << endl;

    slicedSurfaceVectorField Sf
    (
        IOobject
        (
            "Sf2",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimArea,
        mesh.faceAreas()
    );

    //Info<< Sf << endl;

    return 0;
}


// ************************************************************************* //
