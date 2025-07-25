/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
    tetgenToFoam

Group
    grpMeshConversionUtilities

Description
    Convert tetgen .ele and .node and .face files to an OpenFOAM mesh.

    Make sure to use add boundary attributes to the smesh file
    (5 fifth column in the element section)
    and run tetgen with -f option.

    Sample smesh file:
    \verbatim
        # cube.smesh -- A 10x10x10 cube
        8 3
        1       0 0 0
        2       0 10 0
        3       10 10 0
        4       10 0 0
        5       0 0 10
        6       0 10 10
        7       10 10 10
        8       10 0 10
        6 1                 # 1 for boundary info present
        4       1 2 3 4 11  # region number 11
        4       5 6 7 8 21  # region number 21
        4       1 2 6 5 3
        4       4 3 7 8 43
        4       1 5 8 4 5
        4       2 6 7 3 65
        0
        0
    \endverbatim

Note
    - for some reason boundary faces point inwards. I just reverse them
      always. Might use some geometric check instead.
    - marked faces might not actually be boundary faces of mesh.
      This is hopefully handled now by first constructing without boundaries
      and then reconstructing with boundary faces.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "cellModel.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Find label of face.
label findFace(const primitiveMesh& mesh, const face& f)
{
    const labelList& pFaces = mesh.pointFaces()[f[0]];

    forAll(pFaces, i)
    {
        label facei = pFaces[i];

        if (mesh.faces()[facei] == f)
        {
            return facei;
        }
    }

    FatalErrorInFunction
        << "Cannot find face " << f << " in mesh." << abort(FatalError);

    return -1;
}



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Convert tetgen .ele and .node and .face files to an OpenFOAM mesh"
    );

    argList::addArgument("prefix", "The prefix for the input tetgen files");
    argList::addBoolOption
    (
        "noFaceFile",
        "Skip reading .face file for boundary information"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    const auto prefix = args.get<fileName>(1);
    const bool readFaceFile = !args.found("noFaceFile");

    const fileName nodeFile(prefix + ".node");
    const fileName eleFile(prefix + ".ele");
    const fileName faceFile(prefix + ".face");

    if (!readFaceFile)
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << endl;
    }
    else
    {
        Info<< "Files:" << endl
            << "    nodes : " << nodeFile << endl
            << "    elems : " << eleFile << endl
            << "    faces : " << faceFile << endl
            << endl;

        Info<< "Reading .face file for boundary information" << nl << endl;
    }

    if (!isFile(nodeFile) || !isFile(eleFile))
    {
        FatalErrorInFunction
            << "Cannot read " << nodeFile << " or " << eleFile
            << exit(FatalError);
    }

    if (readFaceFile && !isFile(faceFile))
    {
        FatalErrorInFunction
            << "Cannot read " << faceFile << endl
            << "Did you run tetgen with -f option?" << endl
            << "If you don't want to read the .face file and thus not have"
            << " patches please\nrerun with the -noFaceFile option"
            << exit(FatalError);
    }


    IFstream nodeStream(nodeFile);

    //
    // Read nodes.
    //

    // Read header.
    string line;

    do
    {
        nodeStream.getLine(line);
    }
    while (line.starts_with('#'));

    IStringStream nodeLine(line);

    label nNodes, nDims, nNodeAttr;
    bool hasRegion;

    nodeLine >> nNodes >> nDims >> nNodeAttr >> hasRegion;


    Info<< "Read .node header:" << endl
        << "    nodes     : " << nNodes << endl
        << "    nDims     : " << nDims << endl
        << "    nAttr     : " << nNodeAttr << endl
        << "    hasRegion : " << hasRegion << endl
        << endl;

    //
    // read points
    //

    pointField points(nNodes);
    Map<label> nodeToPoint(2*nNodes);

    {
        labelList pointIndex(nNodes);

        label pointi = 0;

        while (nodeStream.good())
        {
            nodeStream.getLine(line);

            if (line.size() && line[0] != '#')
            {
                IStringStream nodeLine(line);

                label nodeI;
                scalar x, y, z;
                label dummy;

                nodeLine >> nodeI >> x >> y >> z;

                for (label i = 0; i < nNodeAttr; i++)
                {
                    nodeLine >> dummy;
                }

                if (hasRegion)
                {
                    nodeLine >> dummy;
                }

                // Store point and node number.
                points[pointi] = point(x, y, z);
                nodeToPoint.insert(nodeI, pointi);
                pointi++;
            }
        }
        if (pointi != nNodes)
        {
            FatalIOErrorInFunction(nodeStream)
                << "Only " << pointi << " nodes present instead of " << nNodes
                << " from header." << exit(FatalIOError);
        }
    }

    //
    // read elements
    //

    IFstream eleStream(eleFile);

    do
    {
        eleStream.getLine(line);
    }
    while (line.starts_with('#'));

    IStringStream eleLine(line);

    label nTets, nPtsPerTet, nElemAttr;

    eleLine >> nTets >> nPtsPerTet >> nElemAttr;


    Info<< "Read .ele header:" << endl
        << "    tets         : " << nTets << endl
        << "    pointsPerTet : " << nPtsPerTet << endl
        << "    nAttr        : " << nElemAttr << endl
        << endl;

    if (nPtsPerTet != 4)
    {
        FatalIOErrorInFunction(eleStream)
            << "Cannot handle tets with "
            << nPtsPerTet << " points per tetrahedron in .ele file" << endl
            << "Can only handle tetrahedra with four points"
            << exit(FatalIOError);
    }

    if (nElemAttr != 0)
    {
        WarningInFunction
            << "Element attributes (third elemenent in .ele header)"
            << " not used" << endl;
    }


    const cellModel& tet = cellModel::ref(cellModel::TET);

    labelList tetPoints(4);

    cellShapeList cells(nTets);
    label celli = 0;

    while (eleStream.good())
    {
        eleStream.getLine(line);

        if (line.size() && line[0] != '#')
        {
            IStringStream eleLine(line);

            label elemI;
            eleLine >> elemI;

            for (label i = 0; i < 4; i++)
            {
                label nodeI;
                eleLine >> nodeI;
                tetPoints[i] = nodeToPoint[nodeI];
            }

            cells[celli++].reset(tet, tetPoints);

            // Skip attributes
            for (label i = 0; i < nElemAttr; i++)
            {
                label dummy;

                eleLine >> dummy;
            }
        }
    }


    //
    // Construct mesh with default boundary only
    //

    auto meshPtr = autoPtr<polyMesh>::New
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.constant(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pointField(points),  // Copy of points
        cells,
        faceListList(),
        wordList(),          // boundaryPatchNames
        wordList(),          // boundaryPatchTypes
        "defaultFaces",
        polyPatch::typeName,
        wordList()
    );
    const polyMesh& mesh = *meshPtr;


    if (readFaceFile)
    {
        label nPatches = 0;

        // List of Foam vertices per boundary face
        faceList boundaryFaces;

        // For each boundary faces the Foam patchID
        labelList boundaryPatch;

        //
        // read boundary faces
        //

        IFstream faceStream(faceFile);

        do
        {
            faceStream.getLine(line);
        }
        while (line.starts_with('#'));

        IStringStream faceLine(line);

        label nFaces, nFaceAttr;

        faceLine >> nFaces >> nFaceAttr;


        Info<< "Read .face header:" << endl
            << "    faces : " << nFaces << endl
            << "    nAttr : " << nFaceAttr << endl
            << endl;


        if (nFaceAttr != 1)
        {
            FatalIOErrorInFunction(faceStream)
                << "Expect boundary markers to be"
                << " present in .face file." << endl
                << "This is the second number in the header which is now:"
                << nFaceAttr << exit(FatalIOError);
        }

        // List of Foam vertices per boundary face
        boundaryFaces.setSize(nFaces);

        // For each boundary faces the Foam patchID
        boundaryPatch.setSize(nFaces);
        boundaryPatch = -1;

        label facei = 0;

        // Region to patch conversion
        Map<label> regionToPatch;

        face f(3);

        while (faceStream.good())
        {
            faceStream.getLine(line);

            if (line.size() && line[0] != '#')
            {
                IStringStream faceLine(line);

                label tetGenFacei, dummy, region;

                faceLine >> tetGenFacei;

                // Read face and reverse orientation (Foam needs outwards
                // pointing)
                for (label i = 0; i < 3; i++)
                {
                    label nodeI;
                    faceLine >> nodeI;
                    f[2-i] = nodeToPoint[nodeI];
                }


                if (findFace(mesh, f) >= mesh.nInternalFaces())
                {
                    boundaryFaces[facei] = f;

                    if (nFaceAttr > 0)
                    {
                        // First attribute is the region number
                        faceLine >> region;


                        // Get Foam patchID and update region->patch table.
                        label patchi = regionToPatch.lookup(region, -1);

                        if (patchi < 0)
                        {
                            patchi = nPatches;
                            regionToPatch.insert(region, nPatches++);

                            Info<< "Mapping tetgen region " << region
                                << " to patch "
                                << patchi << endl;
                        }

                        boundaryPatch[facei] = patchi;

                        // Skip remaining attributes
                        for (label i = 1; i < nFaceAttr; i++)
                        {
                            faceLine >> dummy;
                        }
                    }

                    facei++;
                }
            }
        }


        // Trim
        boundaryFaces.setSize(facei);
        boundaryPatch.setSize(facei);


        // Print region to patch mapping
        Info<< "Regions:" << endl;

        forAllConstIters(regionToPatch, iter)
        {
            Info<< "    region:" << iter.key()
                << '\t' << "patch:" << iter.val() << nl;
        }
        Info<< endl;


        // Storage for boundary faces
        faceListList patchFaces(nPatches);
        wordList patchNames(nPatches);

        forAll(patchNames, patchi)
        {
            patchNames[patchi] = polyPatch::defaultName(patchi);
        }

        wordList patchTypes(nPatches, polyPatch::typeName);
        word defaultFacesName = "defaultFaces";
        word defaultFacesType = polyPatch::typeName;
        wordList patchPhysicalTypes(nPatches, polyPatch::typeName);


        // Sort boundaryFaces by patch using boundaryPatch.
        List<DynamicList<face>> allPatchFaces(nPatches);

        forAll(boundaryPatch, facei)
        {
            label patchi = boundaryPatch[facei];

            allPatchFaces[patchi].append(boundaryFaces[facei]);
        }

        Info<< "Patch sizes:" << endl;

        forAll(allPatchFaces, patchi)
        {
            Info<< "    " << patchNames[patchi] << " : "
                << allPatchFaces[patchi].size() << endl;

            patchFaces[patchi].transfer(allPatchFaces[patchi]);
        }

        Info<< endl;


        meshPtr.reset
        (
            new polyMesh
            (
                IOobject
                (
                    polyMesh::defaultRegion,
                    runTime.constant(),
                    runTime
                ),
                std::move(points),
                cells,
                patchFaces,
                patchNames,
                patchTypes,
                defaultFacesName,
                defaultFacesType,
                patchPhysicalTypes
            )
        );
    }

    // More precision (for points data)
    IOstream::minPrecision(10);

    Info<< "Writing mesh to " << runTime.constant() << endl << endl;

    meshPtr().removeFiles();
    meshPtr().write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
