/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    refineMesh

Group
    grpMeshManipulationUtilities

Description
    Utility to refine cells in multiple directions.

    Command-line option handling:
    - If -all specified or no refineMeshDict exists or, refine all cells
    - If -dict \<file\> specified refine according to \<file\>
    - If refineMeshDict exists refine according to refineMeshDict

    When the refinement or all cells is selected apply 3D refinement for 3D
    cases and 2D refinement for 2D cases.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "polyMesh.H"
#include "Time.H"
#include "cellSet.H"
#include "multiDirRefinement.H"
#include "labelIOList.H"
#include "IOdictionary.H"
#include "syncTools.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Max cos angle for edges to be considered aligned with axis.
static const scalar edgeTol = 1e-3;


// Print edge statistics on mesh.
void printEdgeStats(const polyMesh& mesh)
{
    label nX = 0;
    label nY = 0;
    label nZ = 0;

    scalar minX = GREAT;
    scalar maxX = -GREAT;
    static const vector x(1, 0, 0);

    scalar minY = GREAT;
    scalar maxY = -GREAT;
    static const vector y(0, 1, 0);

    scalar minZ = GREAT;
    scalar maxZ = -GREAT;
    static const vector z(0, 0, 1);

    scalar minOther = GREAT;
    scalar maxOther = -GREAT;

    bitSet isMasterEdge(syncTools::getMasterEdges(mesh));

    const edgeList& edges = mesh.edges();

    forAll(edges, edgeI)
    {
        if (isMasterEdge.test(edgeI))
        {
            const edge& e = edges[edgeI];

            vector eVec(e.vec(mesh.points()));

            scalar eMag = mag(eVec);

            eVec /= eMag;

            if (mag(eVec & x) > 1-edgeTol)
            {
                minX = min(minX, eMag);
                maxX = max(maxX, eMag);
                nX++;
            }
            else if (mag(eVec & y) > 1-edgeTol)
            {
                minY = min(minY, eMag);
                maxY = max(maxY, eMag);
                nY++;
            }
            else if (mag(eVec & z) > 1-edgeTol)
            {
                minZ = min(minZ, eMag);
                maxZ = max(maxZ, eMag);
                nZ++;
            }
            else
            {
                minOther = min(minOther, eMag);
                maxOther = max(maxOther, eMag);
            }
        }
    }

    label nEdges = mesh.nEdges();
    reduce(nEdges, sumOp<label>());
    reduce(nX, sumOp<label>());
    reduce(nY, sumOp<label>());
    reduce(nZ, sumOp<label>());

    reduce(minX, minOp<scalar>());
    reduce(maxX, maxOp<scalar>());

    reduce(minY, minOp<scalar>());
    reduce(maxY, maxOp<scalar>());

    reduce(minZ, minOp<scalar>());
    reduce(maxZ, maxOp<scalar>());

    reduce(minOther, minOp<scalar>());
    reduce(maxOther, maxOp<scalar>());


    Info<< "Mesh edge statistics:" << nl
        << "    x aligned :  number:" << nX << "\tminLen:" << minX
        << "\tmaxLen:" << maxX << nl
        << "    y aligned :  number:" << nY << "\tminLen:" << minY
        << "\tmaxLen:" << maxY << nl
        << "    z aligned :  number:" << nZ << "\tminLen:" << minZ
        << "\tmaxLen:" << maxZ << nl
        << "    other     :  number:" << nEdges - nX - nY - nZ
        << "\tminLen:" << minOther
        << "\tmaxLen:" << maxOther << nl << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Refine cells in multiple directions"
    );

    #include "addOverwriteOption.H"
    #include "addRegionOption.H"

    argList::addOption("dict", "file", "Alternative refineMeshDict");

    argList::addBoolOption
    (
        "all",
        "Refine all cells"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    printEdgeStats(mesh);

    //
    // Read/construct control dictionary
    //

    const bool refineAllCells = args.found("all");
    const bool overwrite = args.found("overwrite");

    // List of cells to refine
    labelList refCells;

    // Dictionary to control refinement
    dictionary refineDict;

    // The -all option has precedence over -dict, or anything else
    if (!refineAllCells)
    {
        const word dictName("refineMeshDict");

        // Obtain dictPath here for messages
        fileName dictPath = args.getOrDefault<fileName>("dict", "");

        IOobject dictIO = IOobject::selectIO
        (
            IOobject
            (
                dictName,
                runTime.system(),
                mesh,
                IOobject::MUST_READ
            ),
            dictPath
        );

        // The reported dictPath will not be full resolved for the output
        // (it will just be the -dict value) but this is purely cosmetic.

        if (dictIO.typeHeaderOk<IOdictionary>(true))
        {
            refineDict = IOdictionary(dictIO);

            Info<< "Refining according to ";

            if (dictPath.empty())
            {
                Info<< dictName;
            }
            else
            {
                Info<< dictPath;
            }
            Info<< nl << endl;
        }
        else if (dictPath.empty())
        {
            Info<< "Refinement dictionary " << dictName << " not found" << nl;
        }
        else
        {
            FatalErrorInFunction
                << "Cannot open specified refinement dictionary "
                << dictPath
                << exit(FatalError);
        }
    }

    if (refineDict.size())
    {
        const word setName(refineDict.get<word>("set"));

        cellSet cells(mesh, setName);

        Info<< "Read " << returnReduce(cells.size(), sumOp<label>())
            << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << endl << endl;

        refCells = cells.toc();
    }
    else
    {
        Info<< "Refining all cells" << nl << endl;

        // Select all cells
        refCells = identity(mesh.nCells());

        if (mesh.nGeometricD() == 3)
        {
            Info<< "3D case; refining all directions" << nl << endl;

            wordList directions(3);
            directions[0] = "tan1";
            directions[1] = "tan2";
            directions[2] = "normal";
            refineDict.add("directions", directions);

            // Use hex cutter
            refineDict.add("useHexTopology", "true");
        }
        else
        {
            const Vector<label> dirs(mesh.geometricD());

            wordList directions(2);

            if (dirs.x() == -1)
            {
                Info<< "2D case; refining in directions y,z\n" << endl;
                directions[0] = "tan2";
                directions[1] = "normal";
            }
            else if (dirs.y() == -1)
            {
                Info<< "2D case; refining in directions x,z\n" << endl;
                directions[0] = "tan1";
                directions[1] = "normal";
            }
            else
            {
                Info<< "2D case; refining in directions x,y\n" << endl;
                directions[0] = "tan1";
                directions[1] = "tan2";
            }

            refineDict.add("directions", directions);

            // Use standard cutter
            refineDict.add("useHexTopology", "false");
        }

        refineDict.add("coordinateSystem", "global");

        dictionary coeffsDict;
        coeffsDict.add("tan1", vector(1, 0, 0));
        coeffsDict.add("tan2", vector(0, 1, 0));
        refineDict.add("globalCoeffs", coeffsDict);

        refineDict.add("geometricCut", "false");
        refineDict.add("writeMesh", "false");
    }


    string oldTimeName(runTime.timeName());

    if (!overwrite)
    {
        ++runTime;
    }


    // Multi-directional refinement (does multiple iterations)
    multiDirRefinement multiRef(mesh, refCells, refineDict);


    // Write resulting mesh
    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }
    mesh.write();


    // Get list of cell splits.
    // (is for every cell in old mesh the cells they have been split into)
    const labelListList& oldToNew = multiRef.addedCells();


    // Create cellSet with added cells for easy inspection
    cellSet newCells(mesh, "refinedCells", refCells.size());

    for (const labelList& added : oldToNew)
    {
        newCells.insert(added);
    }

    Info<< "Writing refined cells ("
        << returnReduce(newCells.size(), sumOp<label>())
        << ") to cellSet "
        << newCells.instance()/newCells.local()/newCells.name()
        << endl << endl;

    newCells.write();


    //
    // Invert cell split to construct map from new to old
    //

    labelIOList newToOld
    (
        IOobject
        (
            "cellMap",
            runTime.timeName(),
            polyMesh::meshSubDir,
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh.nCells()
    );
    newToOld.note() =
        "From cells in mesh at "
      + runTime.timeName()
      + " to cells in mesh at "
      + oldTimeName;


    forAll(oldToNew, oldCelli)
    {
        const labelList& added = oldToNew[oldCelli];

        if (added.size())
        {
            for (const label celli : added)
            {
                newToOld[celli] = oldCelli;
            }
        }
        else
        {
            // Unrefined cell
            newToOld[oldCelli] = oldCelli;
        }
    }

    Info<< "Writing map from new to old cell to "
        << newToOld.objectPath() << nl << endl;

    newToOld.write();

    printEdgeStats(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
