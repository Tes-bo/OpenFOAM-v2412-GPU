/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    Check topology of poly patches used for finite-area.

Input
    mesh (polyMesh)
    meshDefDict (system/faMeshDefintion)

\*---------------------------------------------------------------------------*/

// Preliminary checks
{
    typedef typename uindirectPrimitivePatch::surfaceTopo TopoType;

    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    wordRes polyPatchNames;
    meshDefDict.readIfPresent("polyMeshPatches", polyPatchNames);

    const labelList patchIDs
    (
        pbm.indices(polyPatchNames, true)  // useGroups
    );

    label nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        nFaceLabels += pbm[patchi].size();
    }

    labelList faceLabels(nFaceLabels);

    nFaceLabels = 0;
    for (const label patchi : patchIDs)
    {
        for (const label facei : pbm[patchi].range())
        {
            faceLabels[nFaceLabels] = facei;
            ++nFaceLabels;
        }
    }

    uindirectPrimitivePatch pp
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );

    // Non-manifold
    labelHashSet badEdges(pp.nEdges()/20);

    labelHashSet* pointSetPtr = nullptr;
    labelHashSet* badEdgesPtr = &badEdges;

    bool foundError = false;

    if (returnReduceAnd(pp.empty()))
    {
        // Empty
    }
    else if (UPstream::parRun())
    {
        const labelList meshEdges
        (
            pp.meshEdges(mesh.edges(), mesh.pointEdges())
        );

        // Parallel - use mesh edges
        // - no check for point-pinch
        // - no check for consistent orientation (if that is posible to
        //   check?)

        // Count number of edge/face connections (globally)
        labelList nEdgeConnections(mesh.nEdges(), Zero);

        const labelListList& edgeFaces = pp.edgeFaces();

        forAll(edgeFaces, edgei)
        {
            nEdgeConnections[meshEdges[edgei]] = edgeFaces[edgei].size();
        }

        // Synchronise across coupled edges
        syncTools::syncEdgeList
        (
            mesh,
            nEdgeConnections,
            plusEqOp<label>(),
            label(0)            // null value
        );

        label labelTyp = TopoType::MANIFOLD;
        forAll(meshEdges, edgei)
        {
            const label meshEdgei = meshEdges[edgei];
            const label numNbrs = nEdgeConnections[meshEdgei];
            if (numNbrs == 1)
            {
                //if (pointSetPtr) pointSetPtr->insert(mesh.edges()[meshEdgei]);
                labelTyp = max(labelTyp, TopoType::OPEN);
            }
            else if (numNbrs == 0 || numNbrs > 2)
            {
                if (pointSetPtr) pointSetPtr->insert(mesh.edges()[meshEdgei]);
                if (badEdgesPtr) badEdgesPtr->insert(edgei);
                labelTyp = max(labelTyp, TopoType::ILLEGAL);
            }
        }
        reduce(labelTyp, maxOp<label>());

        foundError = (labelTyp == TopoType::ILLEGAL);
    }
    else
    {
        const TopoType pTyp = pp.surfaceType(badEdgesPtr);

        foundError = (pTyp == TopoType::ILLEGAL);
    }

    if (foundError)
    {
        edgeList dumpEdges(pp.edges(), badEdges.sortedToc());

        vtk::lineWriter writer
        (
            pp.localPoints(),
            dumpEdges,
            fileName
            (
                mesh.time().globalPath()
              / ("faMesh-construct.illegalEdges")
            )
        );

        writer.writeGeometry();

        // CellData
        writer.beginCellData();
        writer.writeProcIDs();

        Info<< "Wrote "
            << returnReduce(dumpEdges.size(), sumOp<label>())
            << " bad edges: " << writer.output().name() << nl;
        writer.close();
    }
}

// ************************************************************************* //
