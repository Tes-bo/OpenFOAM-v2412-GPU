/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "checkTopology.H"
#include "polyMesh.H"
#include "Time.H"
#include "regionSplit.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "IOmanip.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "foamVtkLineWriter.H"
#include "vtkCoordSetWriter.H"
#include "vtkSurfaceWriter.H"
#include "checkTools.H"
#include "treeBoundBox.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// NEWER CODE
// ~~~~~~~~~~
// Return true on error
template<class PatchType>
bool Foam::checkPatch
(
    const bool allGeometry,
    const std::string& name,
    const polyMesh& mesh,
    const PatchType& pp,
    const labelUList& meshEdges,
    labelHashSet* pointSetPtr,
    labelHashSet* badEdgesPtr
)
{
    if (badEdgesPtr)
    {
        badEdgesPtr->clear();
    }

    typedef typename PatchType::surfaceTopo TopoType;

    bool foundError = false;

    const label globalSize = returnReduce(pp.size(), sumOp<label>());

    Info<< "    "
        << setw(20) << name.c_str()
        << setw(9) << globalSize
        << setw(9) << returnReduce(pp.nPoints(), sumOp<label>());

    if (globalSize == 0)
    {
        Info<< setw(34) << "ok (empty)";
    }
    else if (UPstream::parRun())
    {
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

        if (labelTyp == TopoType::MANIFOLD)
        {
            Info<< setw(34) << "ok (closed singly connected)";
        }
        else if (labelTyp == TopoType::OPEN)
        {
            Info<< setw(34) << "ok (non-closed singly connected)";
        }
        else
        {
            Info<< setw(34) << "multiply connected (shared edge)";
        }

        foundError = (labelTyp == TopoType::ILLEGAL);
    }
    else
    {
        const TopoType pTyp = pp.surfaceType(badEdgesPtr);

        if (pTyp == TopoType::MANIFOLD)
        {
            if (pp.checkPointManifold(true, pointSetPtr))
            {
                Info<< setw(34) << "multiply connected (shared point)";
            }
            else
            {
                Info<< setw(34) << "ok (closed singly connected)";
            }

            if (pointSetPtr)
            {
                // Add points on non-manifold edges to make set complete
                pp.checkTopology(false, pointSetPtr);
            }
        }
        else
        {
            if (pointSetPtr)
            {
                pp.checkTopology(false, pointSetPtr);
            }

            if (pTyp == TopoType::OPEN)
            {
                Info<< setw(34) << "ok (non-closed singly connected)";
            }
            else
            {
                Info<< setw(34) << "multiply connected (shared edge)";
            }
        }

        foundError = (pTyp == TopoType::ILLEGAL);
    }

    if (allGeometry)
    {
        boundBox bb(pp.box());
        bb.reduce();

        if (bb.good())
        {
            Info<< ' ' << bb;
        }
    }

    return foundError;
}


// OLDER CODE
// ~~~~~~~~~~
// Return true on error
template<class PatchType>
bool Foam::checkPatch
(
    const bool allGeometry,
    const std::string& name,
    const polyMesh& mesh,
    const PatchType& pp,
    const labelUList& meshFaces,
    const labelUList& meshEdges,
    labelHashSet* pointSetPtr,
    labelHashSet* badEdgesPtr
)
{
    if (badEdgesPtr)
    {
        badEdgesPtr->clear();
    }

    typedef typename PatchType::surfaceTopo TopoType;

    bool foundError = false;

    const label globalSize = returnReduce(pp.size(), sumOp<label>());

    Info<< "    "
        << setw(20) << name.c_str()
        << setw(9) << globalSize
        << setw(9) << returnReduce(pp.nPoints(), sumOp<label>());

    if (globalSize == 0)
    {
        Info<< setw(34) << "ok (empty)";
    }
    else if (UPstream::parRun())
    {
        // Parallel - use mesh edges
        // - no check for point-pinch
        // - no check for consistent orientation (if that is posible to
        //   check?)

        // OLDER CODE
        // ~~~~~~~~~~
        // Synchronise connected faces using global face numbering
        //
        // (see addPatchCellLayer::globalEdgeFaces)
        // From mesh edge to global face labels. Non-empty sublists only for
        // pp edges.
        labelListList globalEdgeFaces(mesh.nEdges());

        const labelListList& edgeFaces = pp.edgeFaces();

        // Global numbering
        const globalIndex globalFaces(mesh.nFaces());

        forAll(edgeFaces, edgei)
        {
            const label meshEdgei = meshEdges[edgei];
            const labelList& eFaces = edgeFaces[edgei];

            // Store face and processor as unique tag.
            labelList& globalEFaces = globalEdgeFaces[meshEdgei];
            globalEFaces.resize(eFaces.size());
            forAll(eFaces, i)
            {
                globalEFaces[i] = globalFaces.toGlobal(meshFaces[eFaces[i]]);
            }
            //Pout<< "At edge:" << meshEdgei
            //    << " ctr:" << mesh.edges()[meshEdgei].centre(mesh.points())
            //    << " have eFaces:" << globalEdgeFaces[meshEdgei]
            //    << endl;
        }

        // Synchronise across coupled edges
        syncTools::syncEdgeList
        (
            mesh,
            globalEdgeFaces,
            ListOps::uniqueEqOp<label>(),
            labelList()         // null value
        );

        label labelTyp = TopoType::MANIFOLD;
        forAll(meshEdges, edgei)
        {
            const label meshEdgei = meshEdges[edgei];
            const labelList& globalEFaces = globalEdgeFaces[meshEdgei];
            const label numNbrs = globalEFaces.size();
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

        if (labelTyp == TopoType::MANIFOLD)
        {
            Info<< setw(34) << "ok (closed singly connected)";
        }
        else if (labelTyp == TopoType::OPEN)
        {
            Info<< setw(34) << "ok (non-closed singly connected)";
        }
        else
        {
            Info<< setw(34) << "multiply connected (shared edge)";
        }

        foundError = (labelTyp == TopoType::ILLEGAL);
    }
    else
    {
        const TopoType pTyp = pp.surfaceType(badEdgesPtr);

        if (pTyp == TopoType::MANIFOLD)
        {
            if (pp.checkPointManifold(true, pointSetPtr))
            {
                Info<< setw(34) << "multiply connected (shared point)";
            }
            else
            {
                Info<< setw(34) << "ok (closed singly connected)";
            }

            if (pointSetPtr)
            {
                // Add points on non-manifold edges to make set complete
                pp.checkTopology(false, pointSetPtr);
            }
        }
        else
        {
            if (pointSetPtr)
            {
                pp.checkTopology(false, pointSetPtr);
            }

            if (pTyp == TopoType::OPEN)
            {
                Info<< setw(34) << "ok (non-closed singly connected)";
            }
            else
            {
                Info<< setw(34) << "multiply connected (shared edge)";
            }
        }

        foundError = (pTyp == TopoType::ILLEGAL);
    }

    if (allGeometry)
    {
        boundBox bb(pp.box());
        bb.reduce();

        if (bb.good())
        {
            Info<< ' ' << bb;
        }
    }

    return foundError;
}


template<class Zone>
Foam::label Foam::checkZones
(
    const polyMesh& mesh,
    const ZoneMesh<Zone, polyMesh>& zones,
    topoSet& set
)
{
    labelList zoneID(set.maxSize(mesh), -1);
    for (const auto& zone : zones)
    {
        for (const label elem : zone)
        {
            if
            (
                zoneID[elem] != -1
             && zoneID[elem] != zone.index()
            )
            {
                set.insert(elem);
            }
            zoneID[elem] = zone.index();
        }
    }

    return returnReduce(set.size(), sumOp<label>());
}


Foam::label Foam::checkTopology
(
    const polyMesh& mesh,
    const bool allTopology,
    const bool allGeometry,
    autoPtr<surfaceWriter>& surfWriter,
    autoPtr<coordSetWriter>& setWriter,
    const bool writeBadEdges
)
{
    label noFailedChecks = 0;

    Info<< "Checking topology..." << endl;

    // Check if the boundary definition is unique
    mesh.boundaryMesh().checkDefinition(true);

    // Check that empty patches cover all sides of the mesh
    {
        label nEmpty = 0;
        forAll(mesh.boundaryMesh(), patchi)
        {
            if (isA<emptyPolyPatch>(mesh.boundaryMesh()[patchi]))
            {
                nEmpty += mesh.boundaryMesh()[patchi].size();
            }
        }
        reduce(nEmpty, sumOp<label>());
        const label nCells = returnReduce(mesh.cells().size(), sumOp<label>());

        // These are actually warnings, not errors.
        if (nCells && (nEmpty % nCells))
        {
            Info<< " ***Total number of faces on empty patches"
                << " is not divisible by the number of cells in the mesh."
                << " Hence this mesh is not 1D or 2D."
                << endl;
        }
    }

    // Check if the boundary processor patches are correct
    mesh.boundaryMesh().checkParallelSync(true);

    // Check names of zones are equal
    mesh.cellZones().checkDefinition(true);
    if (mesh.cellZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }
    mesh.faceZones().checkDefinition(true);
    if (mesh.faceZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }
    mesh.pointZones().checkDefinition(true);
    if (mesh.pointZones().checkParallelSync(true))
    {
        noFailedChecks++;
    }


    {
        cellSet cells(mesh, "illegalCells", mesh.nCells()/100);

        forAll(mesh.cells(), celli)
        {
            const cell& cFaces = mesh.cells()[celli];

            if (cFaces.size() <= 3)
            {
                cells.insert(celli);
            }
            for (const label facei : cFaces)
            {
                if (facei < 0 || facei >= mesh.nFaces())
                {
                    cells.insert(celli);
                    break;
                }
            }
        }
        const label nCells = returnReduce(cells.size(), sumOp<label>());

        if (nCells > 0)
        {
            Info<< "    Illegal cells (less than 4 faces or out of range faces)"
                << " found,  number of cells: " << nCells << endl;
            noFailedChecks++;

            Info<< "  <<Writing " << nCells
                << " illegal cells to set " << cells.name() << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, cells);
            }
        }
        else
        {
            Info<< "    Cell to face addressing OK." << endl;
        }
    }


    {
        pointSet points(mesh, "unusedPoints", mesh.nPoints()/100);
        if (mesh.checkPoints(true, &points))
        {
            noFailedChecks++;

            const label nPoints = returnReduce(points.size(), sumOp<label>());

            Info<< "  <<Writing " << nPoints
                << " unused points to set " << points.name() << endl;
            points.instance() = mesh.pointsInstance();
            points.write();
            if (setWriter && setWriter->enabled())
            {
                mergeAndWrite(*setWriter, points);
            }
        }
    }

    {
        faceSet faces(mesh, "upperTriangularFace", mesh.nFaces()/100);
        if (mesh.checkUpperTriangular(true, &faces))
        {
            noFailedChecks++;
        }

        const label nFaces = returnReduce(faces.size(), sumOp<label>());

        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " unordered faces to set " << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    {
        faceSet faces(mesh, "outOfRangeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceVertices(true, &faces))
        {
            noFailedChecks++;

            const label nFaces = returnReduce(faces.size(), sumOp<label>());

            Info<< "  <<Writing " << nFaces
                << " faces with out-of-range or duplicate vertices to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    if (allTopology)
    {
        cellSet cells(mesh, "zipUpCells", mesh.nCells()/100);
        if (mesh.checkCellsZipUp(true, &cells))
        {
            noFailedChecks++;

            const label nCells = returnReduce(cells.size(), sumOp<label>());

            Info<< "  <<Writing " << nCells
                << " cells with over used edges to set " << cells.name()
                << endl;
            cells.instance() = mesh.pointsInstance();
            cells.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, cells);
            }
        }
    }

    if (allTopology)
    {
        faceSet faces(mesh, "edgeFaces", mesh.nFaces()/100);
        if (mesh.checkFaceFaces(true, &faces))
        {
            noFailedChecks++;
        }

        const label nFaces = returnReduce(faces.size(), sumOp<label>());
        if (nFaces > 0)
        {
            Info<< "  <<Writing " << nFaces
                << " faces with non-standard edge connectivity to set "
                << faces.name() << endl;
            faces.instance() = mesh.pointsInstance();
            faces.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, faces);
            }
        }
    }

    if (allTopology)
    {
        labelList nInternalFaces(mesh.nCells(), Zero);

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            nInternalFaces[mesh.faceOwner()[facei]]++;
            nInternalFaces[mesh.faceNeighbour()[facei]]++;
        }
        const polyBoundaryMesh& patches = mesh.boundaryMesh();
        forAll(patches, patchi)
        {
            if (patches[patchi].coupled())
            {
                const labelUList& owners = patches[patchi].faceCells();

                for (const label facei : owners)
                {
                    nInternalFaces[facei]++;
                }
            }
        }

        cellSet oneCells(mesh, "oneInternalFaceCells", mesh.nCells()/100);
        cellSet twoCells(mesh, "twoInternalFacesCells", mesh.nCells()/100);

        forAll(nInternalFaces, celli)
        {
            if (nInternalFaces[celli] <= 1)
            {
                oneCells.insert(celli);
            }
            else if (nInternalFaces[celli] == 2)
            {
                twoCells.insert(celli);
            }
        }

        const label nOneCells = returnReduce(oneCells.size(), sumOp<label>());

        if (nOneCells > 0)
        {
            Info<< "  <<Writing " << nOneCells
                << " cells with zero or one non-boundary face to set "
                << oneCells.name()
                << endl;
            oneCells.instance() = mesh.pointsInstance();
            oneCells.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, oneCells);
            }
        }

        const label nTwoCells = returnReduce(twoCells.size(), sumOp<label>());

        if (nTwoCells > 0)
        {
            Info<< "  <<Writing " << nTwoCells
                << " cells with two non-boundary faces to set "
                << twoCells.name()
                << endl;
            twoCells.instance() = mesh.pointsInstance();
            twoCells.write();
            if (surfWriter && surfWriter->enabled())
            {
                mergeAndWrite(*surfWriter, twoCells);
            }
        }
    }

    {
        regionSplit rs(mesh);

        if (rs.nRegions() <= 1)
        {
            Info<< "    Number of regions: " << rs.nRegions() << " (OK)."
                << endl;

        }
        else
        {
            Info<< "   *Number of regions: " << rs.nRegions() << endl;

            Info<< "    The mesh has multiple regions which are not connected "
                   "by any face." << endl
                << "  <<Writing region information to "
                << mesh.time().timeName()/"cellToRegion"
                << endl;

            IOListRef<label>
            (
                IOobject
                (
                    "cellToRegion",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                rs
            ).write();


            // Points in multiple regions
            pointSet points
            (
                mesh,
                "multiRegionPoints",
                mesh.nPoints()/1000
            );

            // Is region disconnected
            boolList regionDisconnected(rs.nRegions(), true);
            if (allTopology)
            {
                // -1   : not assigned
                // -2   : multiple regions
                // >= 0 : single region
                labelList pointToRegion(mesh.nPoints(), -1);

                for
                (
                    label facei = mesh.nInternalFaces();
                    facei < mesh.nFaces();
                    ++facei
                )
                {
                    const label regioni = rs[mesh.faceOwner()[facei]];
                    const face& f = mesh.faces()[facei];
                    for (const label verti : f)
                    {
                        label& pRegion = pointToRegion[verti];
                        if (pRegion == -1)
                        {
                            pRegion = regioni;
                        }
                        else if (pRegion == -2)
                        {
                            // Already marked
                            regionDisconnected[regioni] = false;
                        }
                        else if (pRegion != regioni)
                        {
                            // Multiple regions
                            regionDisconnected[regioni] = false;
                            regionDisconnected[pRegion] = false;
                            pRegion = -2;
                            points.insert(verti);
                        }
                    }
                }

                Pstream::listCombineReduce(regionDisconnected, andEqOp<bool>());
            }



            // write cellSet for each region
            PtrList<cellSet> cellRegions(rs.nRegions());
            for (label i = 0; i < rs.nRegions(); i++)
            {
                cellRegions.set
                (
                    i,
                    new cellSet
                    (
                        mesh,
                        "region" + Foam::name(i),
                        mesh.nCells()/100
                    )
                );
            }

            forAll(rs, i)
            {
                cellRegions[rs[i]].insert(i);
            }

            for (label i = 0; i < rs.nRegions(); i++)
            {
                Info<< "  <<Writing region " << i;
                if (allTopology)
                {
                    if (regionDisconnected[i])
                    {
                        Info<< " (fully disconnected)";
                    }
                    else
                    {
                        Info<< " (point connected)";
                    }
                }
                Info<< " with "
                    << returnReduce(cellRegions[i].size(), sumOp<scalar>())
                    << " cells to cellSet " << cellRegions[i].name() << endl;

                cellRegions[i].write();
            }

            const label nPoints = returnReduce(points.size(), sumOp<label>());
            if (nPoints)
            {
                Info<< "  <<Writing " << nPoints
                    << " points that are in multiple regions to set "
                    << points.name() << endl;
                points.write();
                if (setWriter && setWriter->enabled())
                {
                    mergeAndWrite(*setWriter, points);
                }
            }
        }
    }

    // Non-manifold points
    pointSet nonManifoldPoints
    (
        mesh,
        "nonManifoldPoints",
        mesh.nPoints()/1000
    );

    {
        Info<< "\nChecking patch topology for multiply connected"
                << " surfaces..." << endl;

        const polyBoundaryMesh& patches = mesh.boundaryMesh();

        Pout.setf(ios_base::left);

        Info<< "    "
            << setw(20) << "Patch"
            << setw(9) << "Faces"
            << setw(9) << "Points"
            << "  Surface topology";
        if (allGeometry)
        {
            Info<< "  Bounding box";
        }
        Info<< endl;

        for (const polyPatch& pp : patches)
        {
            if (!UPstream::parRun() || !isA<processorPolyPatch>(pp))
            {
                checkPatch
                (
                    allGeometry,
                    pp.name(),
                    mesh,
                    pp,
                    // identity(pp.size(), pp.start()),
                    pp.meshEdges(),
                   &nonManifoldPoints
                );
                Info<< endl;
            }
        }

        // All non-processor boundary patches
        {
            const label nGlobalPatches
            (
                UPstream::parRun()
              ? patches.nNonProcessor()
              : patches.size()
            );

            labelList faceLabels
            (
                identity
                (
                    (
                        patches.range(nGlobalPatches-1).end_value()
                      - patches.start()
                    ),
                    patches.start()
                )
            );

            uindirectPrimitivePatch pp
            (
                UIndirectList<face>(mesh.faces(), faceLabels),
                mesh.points()
            );

            // Non-manifold
            labelHashSet badEdges(pp.nEdges()/20);

            bool hadBadEdges = checkPatch
            (
                allGeometry,
                "\".*\"",
                 mesh,
                pp,
                // faceLabels,
                pp.meshEdges(mesh.edges(), mesh.pointEdges()),
                nullptr,  // No point set
               &badEdges
            );
            Info<< nl << endl;

            if (writeBadEdges && hadBadEdges)
            {
                edgeList dumpEdges(pp.edges(), badEdges.sortedToc());

                vtk::lineWriter writer
                (
                    pp.localPoints(),
                    dumpEdges,
                    fileName
                    (
                        mesh.time().globalPath()
                      / ("checkMesh-illegal-edges")
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
            else if (hadBadEdges)
            {
                Info<< "Detected "
                    << returnReduce(badEdges.size(), sumOp<label>())
                    << " bad edges (possibly relevant for finite-area)" << nl;
            }
        }

        //Info.setf(ios_base::right);
    }

    {
        Info<< "\nChecking faceZone topology for multiply connected"
            << " surfaces..." << endl;

        Pout.setf(ios_base::left);

        const faceZoneMesh& faceZones = mesh.faceZones();

        if (faceZones.size())
        {
            Info<< "    "
                << setw(20) << "FaceZone"
                << setw(9) << "Faces"
                << setw(9) << "Points"
                << setw(34) << "Surface topology";
            if (allGeometry)
            {
                Info<< " Bounding box";
            }
            Info<< endl;

            for (const faceZone& fz : faceZones)
            {
                checkPatch
                (
                    allGeometry,
                    fz.name(),
                    mesh,
                    fz(),           // patch
                    // fz,             // mesh face labels
                    fz.meshEdges(), // mesh edge labels
                   &nonManifoldPoints
                );
                Info<< endl;
            }

            // Check for duplicates
            if (allTopology)
            {
                faceSet mzFaces(mesh, "multiZoneFaces", mesh.nFaces()/100);
                const label nMulti = checkZones(mesh, faceZones, mzFaces);
                if (nMulti)
                {
                    Info<< "  <<Writing " << nMulti
                        << " faces that are in multiple zones"
                        << " to set " << mzFaces.name() << endl;
                    mzFaces.instance() = mesh.pointsInstance();
                    mzFaces.write();
                    if (surfWriter && surfWriter->enabled())
                    {
                        mergeAndWrite(*surfWriter, mzFaces);
                    }
                }
            }
        }
        else
        {
            Info<< "    No faceZones found."<<endl;
        }
    }


    const label nPoints =
        returnReduce(nonManifoldPoints.size(), sumOp<label>());

    if (nPoints)
    {
        Info<< "  <<Writing " << nPoints
            << " conflicting points to set " << nonManifoldPoints.name()
            << endl;
        nonManifoldPoints.instance() = mesh.pointsInstance();
        nonManifoldPoints.write();
        if (setWriter && setWriter->enabled())
        {
            mergeAndWrite(*setWriter, nonManifoldPoints);
        }
    }

    {
        Info<< "\nChecking basic cellZone addressing..." << endl;

        Pout.setf(ios_base::left);

        const cellZoneMesh& cellZones = mesh.cellZones();

        if (cellZones.size())
        {
            Info<< "    "
                << setw(20) << "CellZone"
                << setw(13) << "Cells"
                << setw(13) << "Points"
                << setw(13) << "Volume"
                << "BoundingBox" << endl;

            const cellList& cells = mesh.cells();
            const faceList& faces = mesh.faces();
            const scalarField& cellVolumes = mesh.cellVolumes();

            bitSet isZonePoint(mesh.nPoints());

            for (const cellZone& cZone : cellZones)
            {
                boundBox bb;
                isZonePoint.reset();  // clears all bits (reset count)
                scalar v = 0.0;

                for (const label celli : cZone)
                {
                    v += cellVolumes[celli];
                    for (const label facei : cells[celli])
                    {
                        const face& f = faces[facei];
                        for (const label verti : f)
                        {
                            if (isZonePoint.set(verti))
                            {
                                bb.add(mesh.points()[verti]);
                            }
                         }
                    }
                }

                bb.reduce();  // Global min/max

                Info<< "    "
                    << setw(19) << cZone.name()
                    << ' ' << setw(12)
                    << returnReduce(cZone.size(), sumOp<label>())
                    << ' ' << setw(12)
                    << returnReduce(isZonePoint.count(), sumOp<label>())
                    << ' ' << setw(12)
                    << returnReduce(v, sumOp<scalar>())
                    << ' ' << bb << endl;
            }


            // Check for duplicates
            if (allTopology)
            {
                cellSet mzCells(mesh, "multiZoneCells", mesh.nCells()/100);
                const label nMulti = checkZones(mesh, cellZones, mzCells);
                if (nMulti)
                {
                    Info<< "  <<Writing " << nMulti
                        << " cells that are in multiple zones"
                        << " to set " << mzCells.name() << endl;
                    mzCells.instance() = mesh.pointsInstance();
                    mzCells.write();
                    if (surfWriter && surfWriter->enabled())
                    {
                        mergeAndWrite(*surfWriter, mzCells);
                    }
                }
            }
        }
        else
        {
            Info<< "    No cellZones found."<<endl;
        }
    }


    {
        Info<< "\nChecking basic pointZone addressing..." << endl;

        Pout.setf(ios_base::left);

        const pointZoneMesh& pointZones = mesh.pointZones();

        if (pointZones.size())
        {
            Info<< "    "
                << setw(20) << "PointZone"
                << setw(8) << "Points"
                << "BoundingBox" << nl;

            for (const auto& zone : pointZones)
            {
                boundBox bb
                (
                    mesh.points(),
                    static_cast<const labelUList&>(zone),
                    true  // Reduce (global min/max)
                );

                Info<< "    "
                    << setw(20) << zone.name()
                    << setw(8)
                    << returnReduce(zone.size(), sumOp<label>())
                    << bb << endl;
            }


            // Check for duplicates
            if (allTopology)
            {
                pointSet mzPoints(mesh, "multiZonePoints", mesh.nPoints()/100);
                const label nMulti = checkZones(mesh, pointZones, mzPoints);
                if (nMulti)
                {
                    Info<< "  <<Writing " << nMulti
                        << " points that are in multiple zones"
                        << " to set " << mzPoints.name() << endl;
                    mzPoints.instance() = mesh.pointsInstance();
                    mzPoints.write();
                    if (setWriter && setWriter->enabled())
                    {
                        mergeAndWrite(*setWriter, mzPoints);
                    }
                }
            }
        }
        else
        {
            Info<< "    No pointZones found."<<endl;
        }
    }


    // Force creation of all addressing if requested.
    // Errors will be reported as required
    if (allTopology)
    {
        mesh.cells();
        mesh.faces();
        mesh.edges();
        mesh.points();
        mesh.faceOwner();
        mesh.faceNeighbour();
        mesh.cellCells();
        mesh.edgeCells();
        mesh.pointCells();
        mesh.edgeFaces();
        mesh.pointFaces();
        mesh.cellEdges();
        mesh.faceEdges();
        mesh.pointEdges();
    }

    return noFailedChecks;
}


// ************************************************************************* //
