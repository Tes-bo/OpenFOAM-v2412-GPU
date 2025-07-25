/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 Zeljko Tukovic, FSB Zagreb.
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
    makeWave

Group
    grpMeshManipulationUtilities

Description
    Make wave at top boundary patch.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dynamicFvMesh.H"
#include "dynamicMotionSolverFvMesh.H"
#include "velocityLaplacianFvMotionSolver.H"
#include "polyPatchID.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mesh motion and topological mesh changes utility"
    );

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"

    dynamicMotionSolverFvMesh& dynMsh =
        refCast<dynamicMotionSolverFvMesh>(mesh);

    motionSolver& motion =
        const_cast<motionSolver&>(dynMsh.motion());

    velocityLaplacianFvMotionSolver& vlMotion =
        refCast<velocityLaplacianFvMotionSolver>
        (
            motion
        );

    pointVectorField& pointMotionU = vlMotion.pointMotionU();


    // Set motion boundary condition at top patch

    word topPatchName("top");
    polyPatchID topPatch(topPatchName, mesh.boundaryMesh());
    if (!topPatch.active())
    {
        FatalError
            << "Patch name " << topPatchName << " not found."
            << abort(FatalError);
    }
    label topPatchIndex = topPatch.index();

    const vectorField& topPatchPoints =
        mesh.boundaryMesh()[topPatchIndex].localPoints();

    vectorField topPatchU(topPatchPoints.size(), Zero);

    scalar L = 1.0;
    scalar a0 = 0.01;

    forAll(topPatchPoints, pointI)
    {
        const scalar x = topPatchPoints[pointI].x();
        topPatchU[pointI] =
            vector
            (
                0,
                a0*::cos(M_PI*x/L)/runTime.deltaTValue(),
                0
            );
    }

    fixedValuePointPatchVectorField& topPatchPointMeshU =
        refCast<fixedValuePointPatchVectorField>
        (
            const_cast<pointPatchVectorField&>
            (
                pointMotionU.boundaryField()[topPatchIndex]
            )
        );

    topPatchPointMeshU == topPatchU;

    fileName meshDir = polyMesh::meshSubDir;

    pointIOField points
    (
        IOobject
        (
            "points",
            runTime.findInstance(meshDir, "points"),
            meshDir,
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    points = vlMotion.newPoints();

    Info<< "Writing points into directory " << points.path() << nl << endl;
    points.write();

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
