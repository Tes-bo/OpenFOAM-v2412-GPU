/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    splitMesh

Group
    grpMeshManipulationUtilities

Description
    Splits mesh by making internal faces external. Uses attachDetach.

    Generates a meshModifier of the form:

    Splitter
    {
        type                       attachDetach;
        faceZoneName               membraneFaces;
        masterPatchName            masterPatch;
        slavePatchName             slavePatch;
        triggerTimes               runTime.value();
    }

    so will detach at the current time and split all faces in membraneFaces
    into masterPatch and slavePatch (which have to be present but of 0 size)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "faceSet.H"
#include "attachDetach.H"
#include "attachPolyTopoChanger.H"
#include "regionSide.H"
#include "primitivePatch.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find edge between points v0 and v1.
label findEdge(const primitiveMesh& mesh, const label v0, const label v1)
{
    const labelList& pEdges = mesh.pointEdges()[v0];

    forAll(pEdges, pEdgeI)
    {
        label edgeI = pEdges[pEdgeI];

        const edge& e = mesh.edges()[edgeI];

        if (e.otherVertex(v0) == v1)
        {
            return edgeI;
        }
    }

    FatalErrorInFunction
        << "Cannot find edge between mesh points " << v0 << " and " << v1
        << abort(FatalError);

    return -1;
}


// Checks whether patch present
void checkPatch(const polyBoundaryMesh& bMesh, const word& name)
{
    const label patchi = bMesh.findPatchID(name);

    if (patchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find patch " << name << nl
            << "It should be present but of zero size" << endl
            << "Valid patches are " << bMesh.names()
            << exit(FatalError);
    }

    if (bMesh[patchi].size())
    {
        FatalErrorInFunction
            << "Patch " << name << " is present but non-zero size"
            << exit(FatalError);
    }
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Splits mesh by making internal faces external at defined faceSet"
    );

    argList::noParallel();
    argList::noFunctionObjects();  // Never use function objects

    #include "addOverwriteOption.H"

    argList::addArgument("faceSet", "The faces used for splitting");
    argList::addArgument("master", "The master patch name");
    argList::addArgument("slave", "The slave patch name");

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    const word setName = args[1];
    const word masterPatch = args[2];
    const word slavePatch = args[3];
    const bool overwrite = args.found("overwrite");

    // List of faces to split
    faceSet facesSet(mesh, setName);

    Info<< "Read " << facesSet.size() << " faces to split" << endl << endl;


    // Convert into labelList and check

    labelList faces(facesSet.toc());

    forAll(faces, i)
    {
        if (!mesh.isInternalFace(faces[i]))
        {
            FatalErrorInFunction
            << "Face " << faces[i] << " in faceSet " << setName
            << " is not an internal face."
            << exit(FatalError);
        }
    }


    // Check for empty master and slave patches
    checkPatch(mesh.boundaryMesh(), masterPatch);
    checkPatch(mesh.boundaryMesh(), slavePatch);


    //
    // Find 'side' of all faces on splitregion. Uses regionSide which needs
    // set of edges on side of this region. Use PrimitivePatch to find these.
    //

    // Addressing on faces only in mesh vertices.
    primitiveFacePatch fPatch
    (
        faceList
        (
            UIndirectList<face>
            (
                mesh.faces(),
                faces
            )
        ),
        mesh.points()
    );

    const labelList& meshPoints = fPatch.meshPoints();

    // Mark all fence edges : edges on boundary of fPatch but not on boundary
    // of polyMesh
    labelHashSet fenceEdges(fPatch.size());

    const labelListList& allEdgeFaces = fPatch.edgeFaces();

    forAll(allEdgeFaces, patchEdgeI)
    {
        if (allEdgeFaces[patchEdgeI].size() == 1)
        {
            const edge& e = fPatch.edges()[patchEdgeI];

            label edgeI =
                findEdge
                (
                    mesh,
                    meshPoints[e.start()],
                    meshPoints[e.end()]
                );

            fenceEdges.insert(edgeI);
        }
    }

    // Find sides reachable from 0th face of faceSet
    label startFacei = faces[0];

    regionSide regionInfo
    (
        mesh,
        facesSet,
        fenceEdges,
        mesh.faceOwner()[startFacei],
        startFacei
    );

    // Determine flip state for all faces in faceSet
    boolList zoneFlip(faces.size());

    forAll(faces, i)
    {
        zoneFlip[i] = !regionInfo.sideOwner().found(faces[i]);
    }


    // Create and add face zones and mesh modifiers
    List<pointZone*> pz(0);
    List<faceZone*> fz(1);
    List<cellZone*> cz(0);

    fz[0] =
        new faceZone
        (
            "membraneFaces",
            std::move(faces),
            std::move(zoneFlip),
            0,
            mesh.faceZones()
        );

    Info<< "Adding point and face zones" << endl;
    mesh.addZones(pz, fz, cz);

    attachPolyTopoChanger splitter(mesh);
    splitter.setSize(1);

    // Add the sliding interface mesh modifier to start working at current
    // time
    splitter.set
    (
        0,
        new attachDetach
        (
            "Splitter",
            0,
            splitter,
            "membraneFaces",
            masterPatch,
            slavePatch,
            scalarField(1, runTime.value())
        )
    );

    Info<< nl << "Constructed topologyModifier:" << endl;
    splitter[0].writeDict(Info);

    if (!overwrite)
    {
        ++runTime;
    }

    splitter.attach();

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    else
    {
        mesh.setInstance(runTime.timeName());
    }

    Info<< "Writing mesh to " << runTime.timeName() << endl;
    if (!mesh.write())
    {
        FatalErrorInFunction
            << "Failed writing polyMesh."
            << exit(FatalError);
    }
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);


    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
