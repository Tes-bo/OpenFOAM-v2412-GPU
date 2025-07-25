/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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
    Test-surfaceTree

Description
    Simple tests for building indexedOctree etc.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clockTime.H"
#include "MeshedSurfaces.H"
#include "indexedOctree.H"
#include "AABBTree.H"
#include "treeDataPrimitivePatch.H"
#include "Random.H"
#include "ListListOps.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Test octree building etc"
    );

    argList::noBanner();
    argList::noParallel();
    argList::addArgument("input", "The input surface file");

    argList::addOption("maxLevel", "int", "The max level");
    argList::addOption("leafSize", "int", "The min leaf size");
    argList::addBoolOption("AABB", "AABBTree instead of indexedOctree");

    argList args(argc, argv);

    const auto importName = args.get<fileName>(1);

    const bool useAABB = args.found("AABB");

    label maxLevel = 10;
    label leafSize = 10;

    if (useAABB)
    {
        // May want different settings..
    }

    args.readIfPresent("maxLevel", maxLevel);
    args.readIfPresent("leafSize", leafSize);

    meshedSurface surf(importName);

    Random rndGen(123456);

    treeBoundBox overallBb
    (
        treeBoundBox(surf.box()).extend(rndGen, 1e-4, ROOTVSMALL)
    );

    Info<< "Surface with " << surf.size() << " faces, "
        << surf.nPoints() << " points" << endl;


    clockTime timing;

    if (useAABB)
    {
        AABBTree<face> tree
        (
            surf,
            surf.points(),
            true,       // Equal bins
            maxLevel,   // maxLevel
            leafSize    // minLeafSize
        );

        Info<< "Built AABBTree, maxLevel:" << maxLevel
            << " minLeaf:" << leafSize
            << " with " << tree.boundBoxes().size()
            << " leaves - " << timing.elapsedTime() << 's' << endl;

        {
            OFstream os("AABBTree.obj");
            tree.writeOBJ(os);

            Info<< "Wrote " << os.name() << endl;
        }

        // Info<< "sizes: ";
        // ListListOps::subSizes
        // (
        //     tree.addressing(),
        //     identityOp{}
        // ).writeList(Info) << endl;
    }
    else
    {
        indexedOctree<treeDataPrimitivePatch<meshedSurface>> tree
        (
            treeDataPrimitivePatch<meshedSurface>(surf, 1e-6),
            overallBb,
            maxLevel,   // maxLevel
            leafSize,   // leafSize
            3.0         // duplicity
        );

        Info<< "Built octree, maxLevel:" << maxLevel
            << " minLeaf:" << leafSize
            << " with " << tree.nodes().size()
            << " nodes, " << tree.nLeafs()
            << " leaves - " << timing.elapsedTime() << 's' << endl;

        {
            OFstream os("indexedOctree.obj");
            tree.writeOBJ(os);

            Info<< "Wrote " << os.name() << endl;
        }
    }


    Info<< "\nEnd\n" << endl;
    return 0;
}


// ************************************************************************* //
