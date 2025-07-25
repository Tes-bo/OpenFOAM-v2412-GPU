/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "adiosCoreWrite.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Member Functions * * * * * * * * * * * * * * * //

void Foam::adiosFoam::adiosCoreWrite::putPatchAttributes(const polyMesh& mesh)
{
    if (UPstream::master())
    {
        const fileName varPath(adiosFoam::regionPath(mesh.name()));

        // Info<<"put patch attributes for " << varPath << nl;

        // Attributes are global - only passed via the master
        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        // Only consider 'real' (non-processor) patches
        const label nPatches = patches.nNonProcessor();

        stringList pNames(nPatches);
        stringList pTypes(nPatches);

        for (label patchi=0; patchi < nPatches; ++patchi)
        {
            const polyPatch& p = patches[patchi];
            pNames[patchi] = p.name();
            pTypes[patchi] = p.type();
        }

        putIntAttribute(varPath / "nPatches",     nPatches);
        putListAttribute(varPath / "patch-names", pNames);
        putListAttribute(varPath / "patch-types", pTypes);
    }
}


void Foam::adiosFoam::adiosCoreWrite::writeMeshPoints(const polyMesh& mesh)
{
    const fileName varPath = adiosFoam::meshPath(mesh.name());

    // polyMesh/nPoints (summary)
    // polyMesh/points: 2D array (N points x 3 coordinates)
    putLabelVariable(varPath / "nPoints", mesh.nPoints());
    putListVariable(varPath / "points",   mesh.points());
}


void Foam::adiosFoam::adiosCoreWrite::writeMeshFaces(const polyMesh& mesh)
{
    const fileName varPath = adiosFoam::meshPath(mesh.name());

    // polyMesh/nCells
    // polyMesh/nFaces
    // polyMesh/nInternalFaces

    putLabelVariable(varPath / "nCells", mesh.nCells());
    putLabelVariable(varPath / "nFaces", mesh.nFaces());
    putLabelVariable(varPath / "nInternalFaces", mesh.nInternalFaces());

    {
        const fileName facesPath = varPath/"faces";

        // polyMesh/faces - save in compact form
        // for this we use two separate lists

        // Not particularly elegant, but should be stable

        const faceList& faces = mesh.faces();
        const label nFaces = faces.size();

        labelList start(nFaces+1);

        // Start as per CompactIOList.C:
        start[0] = 0;
        forAll(faces, facei)
        {
            label prev = start[facei];
            start[facei+1] = prev + faces[facei].size();

            if (start[facei+1] < prev)
            {
                FatalErrorInFunction
                    << "Overall number of elements " << start[facei+1]
                    << " of CompactIOList of size " << nFaces
                    << " overflows the representation of a label"
                    << endl
                    << "Please recompile with a larger representation"
                    << " for label" << exit(FatalIOError);
            }
        }

        putListVariable(facesPath/"indices", start);

        labelList elems(start[start.size()-1]);

        // Content as per CompactIOList.C:
        label elemi = 0;
        for (const face& f : faces)
        {
            for (const label pointi : f)
            {
                elems[elemi++] = pointi;
            }
        }

        putListVariable(facesPath/"content", elems);
    }


    // polyMesh/owner - direct write
    putListVariable(varPath/"owner", mesh.faceOwner());

    // polyMesh/neighbour - direct write
    putListVariable(varPath/"neighbour", mesh.faceNeighbour());

    // polyMesh/boundary - byte-stream
    putStreamVariable(varPath/"boundary", mesh.boundaryMesh());
}


void Foam::adiosFoam::adiosCoreWrite::writeMeshAddressing(const polyMesh& mesh)
{
    const fileName varPath = adiosFoam::meshPath(mesh.name());

    // processor addressing types
    //
    // polyMesh/boundaryProcAddressing
    // polyMesh/cellProcAddressing
    // polyMesh/faceProcAddressing
    // polyMesh/pointProcAddressing

    // Q: Do we wish to handle cases with bad addressing?
    // - ie, sanity check that the number of cells/boundaries etc make sense?
    //
    // Might be best to leave this to the reader-side of things instead.

    const fileName addrPath = varPath/"procAddressing";

    const word facesInst = mesh.facesInstance();

    for (const char*  category : { "boundary", "cell", "face", "point" })
    {
        labelIOList list
        (
            IOobject
            (
                word(category) + "ProcAddressing",
                facesInst,
                polyMesh::meshSubDir,
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        );

        // polyMesh/procAddressing/XXX - direct write
        putListVariable(addrPath/category, list);
    }
}


void Foam::adiosFoam::adiosCoreWrite::writeMesh(const polyMesh& mesh)
{
    const fileName varPath = adiosFoam::meshPath(mesh.name());

    writeMeshPoints(mesh);
    writeMeshFaces(mesh);
    // writeMeshAddressing(mesh);

    // zones
    writeMeshZones(varPath/"cellZones",  mesh.cellZones());
    writeMeshZones(varPath/"faceZones",  mesh.faceZones());
    writeMeshZones(varPath/"pointZones", mesh.pointZones());

    // TODO?
    // perhaps not really useful
    // sets/ --> writeMeshSets(mesh);
}


// ************************************************************************* //
