/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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
    refineWallLayer

Group
    grpMeshAdvancedUtilities

Description
    Refine cells next to specified patches.

    Arguments:
        1: List of patch names or regular expressions
        2: The size of the refined cells as a fraction of the edge-length.

    Examples:
        Split the near-wall cells of patch Wall in the middle
            refineWallLayer "(Wall)" 0.5

        Split the near-wall cells of patches Wall1 and Wall2 in the middle
            refineWallLayer "(Wall1 Wall2)" 0.5

        Split the near-wall cells of all patches with names beginning with wall
        with the near-wall cells 10% of the thickness of the original cells
            refineWallLayer '("Wall.*")' 0.1

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "polyTopoChange.H"
#include "cellCuts.H"
#include "cellSet.H"
#include "meshCutter.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Refine cells next to specified patches."
    );

    #include "addOverwriteOption.H"
    argList::addArgument
    (
        "patches",
        "The list of patch names or regex - Eg, '(top \"Wall.\")'"
    );
    argList::addArgument
    (
        "edgeFraction",
        "The size of the refined cells as a fraction of the edge-length"
        " on a (0,1) interval"
    );

    argList::addOption
    (
        "useSet",
        "name",
        "Restrict cells to refine based on specified cellSet name"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createPolyMesh.H"

    const word oldInstance = mesh.pointsInstance();

    // Find set of patches from the list of regular expressions provided
    const wordRes patches(args.getList<wordRe>(1));
    const scalar weight  = args.get<scalar>(2);
    const bool overwrite = args.found("overwrite");

    const labelHashSet patchSet(mesh.boundaryMesh().patchSet(patches));
    if (!patchSet.size())
    {
        FatalErrorInFunction
            << "Cannot find any patches in set " << patches << endl
            << "Valid patches are " << mesh.boundaryMesh().names()
            << exit(FatalError);
    }

    label nPatchFaces = 0;
    label nPatchEdges = 0;

    for (const label patchi : patchSet)
    {
        nPatchFaces += mesh.boundaryMesh()[patchi].size();
        nPatchEdges += mesh.boundaryMesh()[patchi].nEdges();
    }

    // Construct from estimate for the number of cells to refine
    labelHashSet cutCells(4*nPatchFaces);

    // Construct from total patch edges in selected patches
    DynamicList<label> allCutEdges(nPatchEdges);
    DynamicList<scalar> allCutEdgeWeights(nPatchEdges);

    // Find cells to refine
    for (const label patchi : patchSet)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        for (const label meshPointi : meshPoints)
        {
            const labelList& pCells = mesh.pointCells()[meshPointi];

            cutCells.insert(pCells);
        }
    }

    // Edit list of cells to refine according to specified set
    word setName;
    if (args.readIfPresent("useSet", setName))
    {
        Info<< "Subsetting cells to cut based on cellSet"
            << setName << nl << endl;

        cellSet cells(mesh, setName);

        Info<< "Read " << cells.size() << " cells from cellSet "
            << cells.instance()/cells.local()/cells.name()
            << nl << endl;


        cutCells.retain(cells);

        Info<< "Removed from cells to cut all the ones not in set "
            << setName << nl << endl;
    }

    // Mark all mesh points on patch
    bitSet vertOnPatch(mesh.nPoints());

    for (const label patchi : patchSet)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        vertOnPatch.set(meshPoints);
    }

    for (const label patchi : patchSet)
    {
        const polyPatch& pp = mesh.boundaryMesh()[patchi];
        const labelList& meshPoints = pp.meshPoints();

        for (const label meshPointi : meshPoints)
        {
            const labelList& pEdges = mesh.pointEdges()[meshPointi];

            for (const label edgei : pEdges)
            {
                const edge& e = mesh.edges()[edgei];

                label otherPointi = e.otherVertex(meshPointi);

                if (!vertOnPatch.test(otherPointi))
                {
                    allCutEdges.append(edgei);

                    if (e.start() == meshPointi)
                    {
                        allCutEdgeWeights.append(weight);
                    }
                    else
                    {
                        allCutEdgeWeights.append(1 - weight);
                    }
                }
            }
        }
    }

    allCutEdges.shrink();
    allCutEdgeWeights.shrink();

    Info<< "Refining:" << nl
        << "    cells:" << cutCells.size() << nl
        << "    edges:" << allCutEdges.size() << endl;

    // Transfer DynamicLists to straight ones.
    scalarField cutEdgeWeights;
    cutEdgeWeights.transfer(allCutEdgeWeights);
    allCutEdgeWeights.clear();


    // Gets cuts across cells from cuts through edges.
    cellCuts cuts
    (
        mesh,
        cutCells.toc(),     // cells candidate for cutting
        labelList(),        // cut vertices
        allCutEdges,        // cut edges
        cutEdgeWeights      // weight on cut edges
    );

    polyTopoChange meshMod(mesh);

    // Cutting engine
    meshCutter cutter(mesh);

    // Insert mesh refinement into polyTopoChange.
    cutter.setRefinement(cuts, meshMod);

    if (!overwrite)
    {
        ++runTime;
    }

    autoPtr<mapPolyMesh> morphMap = meshMod.changeMesh(mesh, false);

    if (morphMap().hasMotionPoints())
    {
        mesh.movePoints(morphMap().preMotionPoints());
    }

    // Update stored labels on meshCutter.
    cutter.updateMesh(morphMap());

    Info<< "Finished refining" << endl;

    if (overwrite)
    {
        mesh.setInstance(oldInstance);
    }

    // Write resulting mesh
    Info<< "Writing refined mesh to time " << runTime.timeName() << endl;

    mesh.write();
    topoSet::removeFiles(mesh);
    processorMeshes::removeFiles(mesh);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
