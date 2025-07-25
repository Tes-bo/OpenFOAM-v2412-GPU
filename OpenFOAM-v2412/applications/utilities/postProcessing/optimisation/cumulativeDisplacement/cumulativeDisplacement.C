/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    cumulativeDisplacement

Description
    Computes and writes the difference between the mesh points in each time
    instance and the ones in the constant folder. Additionally, the projection
    of this difference to the normal point vectors of the initial mesh is also
    written

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "emptyFvPatch.H"
#include "coupledFvPatch.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "addRegionOption.H"

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    #include "createNamedMesh.H"
    #include "createFields.H"

    // polyPatch::pointNormals will give the wrong result for points
    // belonging to multiple patches or patch-processorPatch intersections.
    // Keeping a mesh-wide field to allow easy reduction using syncTools.
    // A bit expensive? Better way?
    vectorField pointNormals(mesh.nPoints(), vector::zero);
    for (const fvPatch& patch : mesh.boundary())
    {
        //  Each local patch point belongs to these local patch faces.
        //  Local numbering
        const labelListList& patchPointFaces = patch.patch().pointFaces();
        if (!isA<coupledFvPatch>(patch) && !isA<emptyFvPatch>(patch))
        {
            const labelList&  meshPoints = patch.patch().meshPoints();
            const vectorField nf(patch.nf());
            forAll(meshPoints, ppI)
            {
                const labelList& pointFaces = patchPointFaces[ppI];
                forAll(pointFaces, pfI)
                {
                    const label& localFaceIndex = pointFaces[pfI];
                    pointNormals[meshPoints[ppI]] += nf[localFaceIndex];
                }
            }
        }
    }
    // Sum from all processors
    syncTools::syncPointList
    (
        mesh, pointNormals, plusEqOp<vector>(), vector::zero
    );
    pointNormals /= mag(pointNormals) + VSMALL;

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        const pointMesh& pMesh(pointMesh::New(mesh));
        // Point displacement projected to the
        // unit point normal of the initial geometry
        pointScalarField normalDisplacement
        (
            IOobject
            (
                "normalDisplacement",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedScalar(dimless, Zero)
        );
        // Point displacement as a vector
        pointVectorField displacement
        (
            IOobject
            (
                "displacement",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            pMesh,
            dimensionedVector(dimless, Zero)
        );

        Info<< "Calculating cumulative mesh movement for time "
            << runTime.timeName() << endl;
        // Normal displacement
        const pointField& meshPoints = mesh.points();
        forAll(mesh.boundary(), pI)
        {
            const polyPatch& patch = mesh.boundaryMesh()[pI];
            const labelList& localPoints = patch.meshPoints();
            forAll(localPoints, ppI)
            {
                label pointI = localPoints[ppI];
                normalDisplacement[pointI]  =
                    (meshPoints[pointI] - points0[pointI])
                  & pointNormals[pointI];
            }
        }
        normalDisplacement.write();
        // Vectorial displacement
        displacement.primitiveFieldRef() = meshPoints - points0;
        displacement.write();
    }

    // Print execution time
    mesh.time().printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
