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
    makePerturbation

Group
    grpMeshManipulationUtilities

Description
    Make droplet perturbation.

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

    // Set motion boundary condition at free-surface patch

    word fsPatchName("freeSurface");
    polyPatchID fsPatch(fsPatchName, mesh.boundaryMesh());
    if (!fsPatch.active())
    {
        FatalError
            << "Patch name " << fsPatchName << " not found."
            << abort(FatalError);
    }
    label fsPatchIndex = fsPatch.index();

    const vectorField& fsPatchPoints =
        mesh.boundaryMesh()[fsPatchIndex].localPoints();

    vectorField fsPatchU(fsPatchPoints.size(), vector::zero);

    // Perturbation amplitude
    scalar epsilon = 0.02;

    // Oscilation mode
    label n = 4;

    forAll(fsPatchPoints, pointI)
    {
        const scalar x = fsPatchPoints[pointI].x();
        const scalar y = fsPatchPoints[pointI].y();

        const scalar theta = ::atan2(y, x);

        scalar scale = (1.0 + epsilon*::cos(n*theta));

        vector r0(x, y, 0);

        fsPatchU[pointI] = (scale - 1.0)*r0/runTime.deltaTValue();
    }

    fixedValuePointPatchVectorField& fsPatchPointMeshU =
        refCast<fixedValuePointPatchVectorField>
        (
            const_cast<pointPatchVectorField&>
            (
                pointMotionU.boundaryField()[fsPatchIndex]
            )
        );

    fsPatchPointMeshU == fsPatchU;

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
