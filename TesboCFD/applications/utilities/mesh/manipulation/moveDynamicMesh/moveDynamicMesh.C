/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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
    moveDynamicMesh

Group
    grpMeshManipulationUtilities

Description
    Mesh motion and topological mesh changes utility.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "dynamicFvMesh.H"
#include "pimpleControl.H"
#include "cyclicAMIPolyPatch.H"
#include "PatchTools.H"
#include "foamVtkSurfaceWriter.H"
#include "functionObject.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Dump patch + weights to vtk file
void writeWeights
(
    const polyMesh& mesh,
    const scalarField& wghtSum,
    const primitivePatch& patch,
    const fileName& directory,
    const fileName& prefix,
    const Time& runTime
)
{
    // Collect geometry
    labelList pointToGlobal;
    labelList uniqueMeshPointLabels;
    autoPtr<globalIndex> globalPoints;
    autoPtr<globalIndex> globalFaces;
    faceList mergedFaces;
    pointField mergedPoints;
    Foam::PatchTools::gatherAndMerge
    (
        mesh,
        patch.localFaces(),
        patch.meshPoints(),
        patch.meshPointMap(),

        pointToGlobal,
        uniqueMeshPointLabels,
        globalPoints,
        globalFaces,

        mergedFaces,
        mergedPoints
    );

    // Collect field
    scalarField mergedWeights;
    globalFaces().gather(wghtSum, mergedWeights);

    instant inst(runTime.value(), runTime.timeName());

    if (Pstream::master())
    {
        vtk::surfaceWriter writer
        (
            mergedPoints,
            mergedFaces,
            (directory/prefix + "_" + inst.name()),
            false // serial: master-only
        );

        writer.setTime(inst);
        writer.writeTimeValue();
        writer.writeGeometry();

        writer.beginCellData(1);
        writer.write("weightsSum", mergedWeights);
    }
}


void writeWeights(const polyMesh& mesh)
{
    const fileName outputDir
    (
        mesh.time().globalPath()/functionObject::outputPrefix/"checkAMI"
    );

    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        const auto* cpp = isA<cyclicAMIPolyPatch>(pp);

        if (cpp && cpp->owner())
        {
            const auto& cycPatch = *cpp;
            const auto& nbrPatch = cycPatch.neighbPatch();

            const AMIPatchToPatchInterpolation& ami = cycPatch.AMI();

            Info<< "Calculating AMI weights between owner patch: "
                << cycPatch.name() << " and neighbour patch: "
                << nbrPatch.name() << endl;

            writeWeights
            (
                mesh,
                ami.tgtWeightsSum(),
                nbrPatch,
                outputDir,
                "patch" + Foam::name(pp.index()) + "-tgt",
                mesh.time()
            );
            writeWeights
            (
                mesh,
                ami.srcWeightsSum(),
                cycPatch,
                outputDir,
                "patch" + Foam::name(pp.index()) + "-src",
                mesh.time()
            );
        }
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Mesh motion and topological mesh changes utility"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"
    argList::addBoolOption
    (
        "checkAMI",
        "Check AMI weights and write VTK files of the AMI patches"
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedDynamicFvMesh.H"

    const bool checkAMI = args.found("checkAMI");

    if (checkAMI)
    {
        Info<< "Writing VTK files with weights of AMI patches." << nl << endl;
    }

    const bool overwrite = args.found("overwrite");
    const word oldInstance = mesh.pointsInstance();


    pimpleControl pimple(mesh);

    bool moveMeshOuterCorrectors
    (
        pimple.dict().getOrDefault("moveMeshOuterCorrectors", false)
    );

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << endl;

        while (pimple.loop())
        {
            if (pimple.firstIter() || moveMeshOuterCorrectors)
            {
                mesh.update();
            }
        }

        if (overwrite)
        {
            mesh.setInstance(oldInstance);
            runTime.write();
            runTime.printExecutionTime(Info);
            break;
        }


        mesh.checkMesh(true);

        if (checkAMI)
        {
            writeWeights(mesh);
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
