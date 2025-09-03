/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "checkTools.H"
#include "polyMesh.H"
#include "globalMeshData.H"
#include "hexMatcher.H"
#include "wedgeMatcher.H"
#include "prismMatcher.H"
#include "pyrMatcher.H"
#include "tetWedgeMatcher.H"
#include "tetMatcher.H"
#include "IOmanip.H"
#include "OFstream.H"
#include "pointSet.H"
#include "faceSet.H"
#include "cellSet.H"
#include "Time.H"
#include "coordSetWriter.H"
#include "surfaceWriter.H"
#include "syncTools.H"
#include "globalIndex.H"
#include "PatchTools.H"
#include "functionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::printMeshStats(const polyMesh& mesh, const bool allTopology)
{
    Info<< "Mesh stats " << mesh.regionName() << nl
        << "    points:           "
        << returnReduce(mesh.points().size(), sumOp<label>()) << nl;


    // Count number of internal points (-1 if not sorted; 0 if no internal
    // points)
    const label minInt = returnReduce(mesh.nInternalPoints(), minOp<label>());
    const label maxInt = returnReduce(mesh.nInternalPoints(), maxOp<label>());

    if (minInt == -1 && maxInt > 0)
    {
        WarningInFunction
            << "Some processors have their points sorted into internal"
            << " and external and some do not." << endl
            << "    This can cause problems later on." << endl;
    }
    else if (minInt != -1)
    {
        // Assume all sorted
        label nInternalPoints = returnReduce
        (
            mesh.nInternalPoints(),
            sumOp<label>()
        );
        Info<< "    internal points:  " << nInternalPoints << nl;
    }

    if (allTopology && (minInt != -1))
    {
        label nEdges = returnReduce(mesh.nEdges(), sumOp<label>());
        label nInternalEdges = returnReduce
        (
            mesh.nInternalEdges(),
            sumOp<label>()
        );
        label nInternal1Edges = returnReduce
        (
            mesh.nInternal1Edges(),
            sumOp<label>()
        );
        label nInternal0Edges = returnReduce
        (
            mesh.nInternal0Edges(),
            sumOp<label>()
        );

        Info<< "    edges:            " << nEdges << nl
            << "    internal edges:   " << nInternalEdges << nl
            << "    internal edges using one boundary point:   "
            << nInternal1Edges-nInternal0Edges << nl
            << "    internal edges using two boundary points:  "
            << nInternalEdges-nInternal1Edges << nl;
    }

    label nFaces = returnReduce(mesh.faces().size(), sumOp<label>());
    label nIntFaces = returnReduce(mesh.faceNeighbour().size(), sumOp<label>());
    label nCells = returnReduce(mesh.cells().size(), sumOp<label>());
    label nPatches = mesh.boundaryMesh().size();

    Info<< "    faces:            " << nFaces << nl
        << "    internal faces:   " << nIntFaces << nl
        << "    cells:            " << nCells << nl
        << "    faces per cell:   "
        << (scalar(nFaces) + scalar(nIntFaces))/max(1, nCells) << nl
        << "    boundary patches: ";

    if (Pstream::parRun())
    {
        // Number of global patches and min-max range of total patches
        Info<< mesh.boundaryMesh().nNonProcessor() << ' '
            << returnReduce(labelMinMax(nPatches), minMaxOp<label>()) << nl;
    }
    else
    {
        Info<< nPatches << nl;
    }

    Info<< "    point zones:      " << mesh.pointZones().size() << nl
        << "    face zones:       " << mesh.faceZones().size() << nl
        << "    cell zones:       " << mesh.cellZones().size() << nl
        << endl;

    // Construct shape recognizers
    prismMatcher prism;
    wedgeMatcher wedge;
    tetWedgeMatcher tetWedge;

    // Counters for different cell types
    label nHex = 0;
    label nWedge = 0;
    label nPrism = 0;
    label nPyr = 0;
    label nTet = 0;
    label nTetWedge = 0;
    label nUnknown = 0;

    Map<label> polyhedralFaces;

    for (label celli = 0; celli < mesh.nCells(); celli++)
    {
        if (hexMatcher::test(mesh, celli))
        {
            nHex++;
        }
        else if (tetMatcher::test(mesh, celli))
        {
            nTet++;
        }
        else if (pyrMatcher::test(mesh, celli))
        {
            nPyr++;
        }
        else if (prism.isA(mesh, celli))
        {
            nPrism++;
        }
        else if (wedge.isA(mesh, celli))
        {
            nWedge++;
        }
        else if (tetWedge.isA(mesh, celli))
        {
            nTetWedge++;
        }
        else
        {
            nUnknown++;
            polyhedralFaces(mesh.cells()[celli].size())++;
        }
    }

    reduce(nHex,sumOp<label>());
    reduce(nPrism,sumOp<label>());
    reduce(nWedge,sumOp<label>());
    reduce(nPyr,sumOp<label>());
    reduce(nTetWedge,sumOp<label>());
    reduce(nTet,sumOp<label>());
    reduce(nUnknown,sumOp<label>());

    Info<< "Overall number of cells of each type:" << nl
        << "    hexahedra:     " << nHex << nl
        << "    prisms:        " << nPrism << nl
        << "    wedges:        " << nWedge << nl
        << "    pyramids:      " << nPyr << nl
        << "    tet wedges:    " << nTetWedge << nl
        << "    tetrahedra:    " << nTet << nl
        << "    polyhedra:     " << nUnknown
        << endl;

    if (nUnknown > 0)
    {
        Pstream::mapCombineGather(polyhedralFaces, plusEqOp<label>());

        Info<< "    Breakdown of polyhedra by number of faces:" << nl
            << "        faces" << "   number of cells" << endl;

        const labelList sortedKeys = polyhedralFaces.sortedToc();

        forAll(sortedKeys, keyi)
        {
            const label nFaces = sortedKeys[keyi];

            Info<< setf(std::ios::right) << setw(13)
                << nFaces << "   " << polyhedralFaces[nFaces] << nl;
        }
    }

    Info<< endl;
}


void Foam::mergeAndWrite
(
    const polyMesh& mesh,
    surfaceWriter& writer,
    const word& name,
    const indirectPrimitivePatch& setPatch,
    const fileName& outputDir
)
{
    writer.open
    (
        setPatch.localPoints(),
        setPatch.localFaces(),
        (outputDir / name)
    );

    writer.write();
    writer.clear();
}


void Foam::mergeAndWrite
(
    surfaceWriter& writer,
    const faceSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), set.sortedToc()),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().globalPath()
      / functionObject::outputPrefix
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();  // Remove unneeded ".."

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    surfaceWriter& writer,
    const cellSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();


    // Determine faces on outside of cellSet
    bitSet isInSet(mesh.nCells());
    for (const label celli : set)
    {
        isInSet.set(celli);
    }


    boolList bndInSet(mesh.nBoundaryFaces());
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        forAll(fc, i)
        {
            bndInSet[pp.start()+i-mesh.nInternalFaces()] = isInSet[fc[i]];
        }
    }
    syncTools::swapBoundaryFaceList(mesh, bndInSet);


    DynamicList<label> outsideFaces(3*set.size());
    for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
    {
        const bool ownVal = isInSet[mesh.faceOwner()[facei]];
        const bool neiVal = isInSet[mesh.faceNeighbour()[facei]];

        if (ownVal != neiVal)
        {
            outsideFaces.append(facei);
        }
    }


    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];
        const labelList& fc = pp.faceCells();
        if (pp.coupled())
        {
            forAll(fc, i)
            {
                label facei = pp.start()+i;

                const bool neiVal = bndInSet[facei-mesh.nInternalFaces()];
                if (isInSet[fc[i]] && !neiVal)
                {
                    outsideFaces.append(facei);
                }
            }
        }
        else
        {
            forAll(fc, i)
            {
                if (isInSet[fc[i]])
                {
                    outsideFaces.append(pp.start()+i);
                }
            }
        }
    }


    const indirectPrimitivePatch setPatch
    (
        IndirectList<face>(mesh.faces(), outsideFaces),
        mesh.points()
    );

    fileName outputDir
    (
        set.time().globalPath()
      / functionObject::outputPrefix
      / mesh.pointsInstance()
      / set.name()
    );
    outputDir.clean();  // Remove unneeded ".."

    mergeAndWrite(mesh, writer, set.name(), setPatch, outputDir);
}


void Foam::mergeAndWrite
(
    coordSetWriter& writer,
    const pointSet& set
)
{
    const polyMesh& mesh = refCast<const polyMesh>(set.db());

    labelField mergedIDs(set.sortedToc());
    pointField mergedPts(mesh.points(), mergedIDs);

    if (Pstream::parRun())
    {
        // Note: we explicitly do not merge the points
        // (mesh.globalData().mergePoints etc) since this might
        // hide any synchronisation problem

        // Renumber local ids -> global ids
        globalIndex(mesh.nPoints()).inplaceToGlobal(mergedIDs);

        globalIndex gatherer(globalIndex::gatherOnly{}, mergedIDs.size());
        gatherer.gatherInplace(mergedIDs);
        gatherer.gatherInplace(mergedPts);
    }


    // Write with pointID
    if (Pstream::master())
    {
        coordSet coords(set.name(), "distance", mergedPts, mag(mergedPts));

        // Output. E.g. pointSet p0 -> postProcessing/<time>/p0.vtk

        fileName outputPath
        (
            set.time().globalPath()
          / functionObject::outputPrefix
          / mesh.pointsInstance()
          / set.name()
        );
        outputPath.clean();  // Remove unneeded ".."

        writer.open(coords, outputPath);
        writer.nFields(1);
        writer.write("pointID", mergedIDs);
        writer.close(true);
    }
}


// ************************************************************************* //
