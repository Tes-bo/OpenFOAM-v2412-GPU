/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM, distributed under GPL-3.0-or-later.

Description
    VTK output of blockMesh topology blocks

    Always write in ASCII since the mesh is small and we want easily
    readable output for inspection.

\*---------------------------------------------------------------------------*/

// Common

refPtr<polyMesh> topoMeshPtr(blocks.topology(true));
const polyMesh& topoMesh = topoMeshPtr();

// Internal mesh - ie, the blocks
{
    const vtk::vtuCells topoCells(topoMesh, vtk::formatType::INLINE_ASCII);

    vtk::internalMeshWriter writer
    (
        topoMesh,
        topoCells,
        vtk::formatType::INLINE_ASCII,
        runTime.path()/"blockTopology"
    );

    Info<< "Writing block topology in vtk format: "
        << args.relativePath(writer.output()).c_str() << endl;

    writer.writeGeometry();
    writer.beginCellData();
    writer.writeCellIDs();

    // No cell decomposition, so there is a 1-to-1 correspondence between
    // input blocks and VTK output cells.

    vectorField localNormal(blocks.size());

    // Generate local normals as fields for visualisation
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        const label faceMin = label(2*cmpt);
        const label faceMax = faceMin+1;

        localNormal.resize(blocks.size());

        forAll(blocks, blocki)
        {
            const cellShape& shape = blocks[blocki].blockShape();
            const pointField& verts = blocks[blocki].vertices();

            if (shape.model() == cellModel::ref(cellModel::HEX))
            {
                localNormal[blocki] =
                (
                    // 30% cell-width as arbitrary value for vector length
                    0.3*mag
                    (
                        shape.face(faceMax).centre(verts)
                      - shape.face(faceMin).centre(verts)
                    )
                  * normalised
                    (
                        // Weigh by area to avoid influence of zero-faces
                        shape.face(faceMax).areaNormal(verts)
                      - shape.face(faceMin).areaNormal(verts)
                    )
                );
            }
            else
            {
                // Could handle other shapes (how?)
                localNormal[blocki] = Zero;
            }
        }

        // Additional safety (should not happen)
        localNormal.resize(topoMesh.nCells(), Zero);

        writer.writeCellData
        (
            word("local-direction" + name(cmpt)),
            localNormal
        );
    }

    // if (topoMesh.nCells() != blocks.size())
    // {
    //     Info<< "Warning: indicated normals may be incorrect" << nl;
    // }
}


// Block boundary faces
{
    const label nIntFaces = topoMesh.nInternalFaces();
    const label nBndFaces = topoMesh.nBoundaryFaces();

    faceList bndFaces
    (
        faceList::subList(topoMesh.faces(), nBndFaces, nIntFaces)
    );

    vtk::surfaceWriter writer
    (
        topoMesh.points(),
        bndFaces,
        vtk::formatType::INLINE_ASCII,
        runTime.path()/"blockFaces"
    );


    labelList blockIds(nBndFaces, -1);
    labelList cellFaceIds(nBndFaces, -1);
    labelList patchIds(nBndFaces, -1);

    {
        const labelList& own = topoMesh.faceOwner();
        const cellList& cells = topoMesh.cells();
        const polyPatchList& patches = topoMesh.boundaryMesh();

        for (const polyPatch& pp : patches)
        {
            label bndFacei = pp.start() - nIntFaces;
            label meshFacei = pp.start();

            forAll(pp, bfacei)
            {
                const label celli = own[meshFacei];
                const label cellFacei = cells[celli].find(meshFacei);

                blockIds[bndFacei] = celli;
                cellFaceIds[bndFacei] = cellFacei;
                patchIds[bndFacei] = pp.index();

                ++bndFacei;
                ++meshFacei;
            }
        }
    }

    Info<< "Writing block boundary faces in vtk format: "
        << args.relativePath(writer.output()).c_str() << endl;

    writer.writeGeometry();

    writer.beginCellData();
    writer.writeCellData("block", blockIds);
    writer.writeCellData("face", cellFaceIds);
    writer.writeCellData("patch", patchIds);
}


// ************************************************************************* //
