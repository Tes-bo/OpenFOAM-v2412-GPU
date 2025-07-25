/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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
    renumberMesh

Group
    grpMeshManipulationUtilities

Description
    Renumbers the cell list in order to reduce the bandwidth, reading and
    renumbering all fields from all the time directories.

    By default uses bandCompression (Cuthill-McKee) or the method
    specified by the -renumber-method option, but will read
    system/renumberMeshDict if -dict option is present

Usage
    \b renumberMesh [OPTIONS]

    Options:
      - \par -allRegions
        Use all regions in regionProperties

      - \par -case \<dir\>
        Specify case directory to use (instead of the cwd).

      - \par -constant
        Include the 'constant/' dir in the times list.

      - \par -decompose
        Aggregate initially with a decomposition method (serial only)

      - \par -decomposeParDict \<file\>
        Use specified file for decomposePar dictionary.

      - \par -dict \<file\>
        Use specified file for renumberMeshDict dictionary.

      - \par -dry-run
        Test only

      - \par -frontWidth
        Calculate the rms of the front-width

      - \par -latestTime
        Select the latest time.

      - \par -lib \<name\>
        Additional library or library list to load (can be used multiple times).

      - \par -no-fields
        Suppress renumber of fields

      - \par -noZero
        Exclude the \a 0 dir from the times list.

      - \par -overwrite
        Overwrite existing mesh/results files

      - \par -parallel
        Run in parallel

      - \par -region \<regionName\>
        Renumber named region.

      - \par -regions \<wordRes\>
        Renumber named regions.

      - \par -renumber-coeffs \<string-content\>
        String to create renumber dictionary contents.

      - \par -renumber-method \<name\>
        Specify renumber method (default: CuthillMcKee) without dictionary

      - \par -time \<value\>
        Specify time to select

      - \par -verbose
        Additional verbosity.

      - \par -doc
        Display documentation in browser.

      - \par -doc-source
        Display source code in browser.

      - \par -help
        Display short help and exit.

      - \par -help-man
        Display full help (manpage format) and exit.

      - \par -help-notes
        Display help notes (description) and exit.

      - \par -help-full
        Display full help and exit.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "clockTime.H"
#include "timeSelector.H"
#include "IOobjectList.H"
#include "fvMesh.H"
#include "polyTopoChange.H"
#include "ReadFields.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SortableList.H"
#include "decompositionMethod.H"
#include "decompositionModel.H"
#include "renumberMethod.H"
#include "foamVtkInternalMeshWriter.H"
#include "CuthillMcKeeRenumber.H"
#include "fvMeshSubset.H"
#include "cellSet.H"
#include "faceSet.H"
#include "pointSet.H"
#include "processorMeshes.H"
#include "hexRef8Data.H"
#include "regionProperties.H"
#include "polyMeshTools.H"
#include "subsetAdjacency.H"

using namespace Foam;

// Slightly messy way of handling timing, but since the timing points
// are scattered between 'main()' and other local functions...

clockTime timer;

// Timing categories
enum TimingType
{
    READ_MESH,    // Reading mesh
    READ_FIELDS,  // Reading fields
    DECOMPOSE,    // Domain decomposition (if any)
    CELL_CELLS,   // globalMeshData::calcCellCells
    RENUMBER,     // The renumberMethod
    REORDER,      // Mesh reordering (topoChange)
    WRITING,      // Writing mesh/fields
};
FixedList<double, 8> timings;


// Create named field from labelList for post-processing
tmp<volScalarField> createScalarField
(
    const fvMesh& mesh,
    const word& name,
    const labelUList& elems
)
{
    auto tfld = volScalarField::New
    (
        name,
        mesh,
        dimensionedScalar(word::null, dimless, -1),
        fvPatchFieldBase::zeroGradientType()
    );
    auto& fld = tfld.ref();

    forAll(elems, celli)
    {
       fld[celli] = elems[celli];
    }

    fld.correctBoundaryConditions();
    return tfld;
}


// Calculate band of mesh
// label getBand(const labelUList& owner, const labelUList& neighbour)
// {
//     label bandwidth = 0;
//
//     forAll(neighbour, facei)
//     {
//         const label width = neighbour[facei] - owner[facei];
//
//         if (bandwidth < width)
//         {
//             bandwidth = width;
//         }
//     }
//     return bandwidth;
// }


// Calculate band and profile of matrix. Profile is scalar to avoid overflow
Tuple2<label, scalar> getBand
(
    const CompactListList<label>& mat
)
{
    Tuple2<label, scalar> metrics(0, 0);

    auto& bandwidth = metrics.first();
    auto& profile = metrics.second();

    forAll(mat, celli)
    {
        const auto& neighbours = mat[celli];

        const label nNbr = neighbours.size();

        if (nNbr)
        {
            // Max distance
            const label width = (neighbours[nNbr-1] - celli);

            if (bandwidth < width)
            {
                bandwidth = width;
            }

            profile += scalar(width);
        }
    }

    return metrics;
}


// Calculate band of matrix
void getBand
(
    const bool calculateIntersect,
    const label nCells,
    const labelUList& owner,
    const labelUList& neighbour,
    label& bandwidth,
    scalar& profile,            // scalar to avoid overflow
    scalar& sumSqrIntersect     // scalar to avoid overflow
)
{
    labelList cellBandwidth(nCells, Foam::zero{});

    bandwidth = 0;

    forAll(neighbour, facei)
    {
        const label own = owner[facei];
        const label nei = neighbour[facei];

        // Note: mag not necessary for correct (upper-triangular) ordering.
        const label width = nei - own;

        if (cellBandwidth[nei] < width)
        {
            cellBandwidth[nei] = width;

            if (bandwidth < width)
            {
                bandwidth = width;
            }
        }
    }

    // Do not use field algebra because of conversion label to scalar
    profile = 0;
    for (const label width : cellBandwidth)
    {
        profile += scalar(width);
    }

    sumSqrIntersect = 0;
    if (calculateIntersect)
    {
        scalarField nIntersect(nCells, Foam::zero{});

        forAll(nIntersect, celli)
        {
            for (label colI = celli-cellBandwidth[celli]; colI <= celli; colI++)
            {
                nIntersect[colI] += scalar(1);
            }
        }

        sumSqrIntersect = sum(Foam::sqr(nIntersect));
    }
}


// Determine upper-triangular face order
labelList getFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder      // New to old cell
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    DynamicList<label> nbr(64);
    DynamicList<label> order(64);

    forAll(cellOrder, newCelli)
    {
        const label oldCelli = cellOrder[newCelli];

        const cell& cFaces = mesh.cells()[oldCelli];

        // Neighbouring cells
        nbr.clear();

        for (const label facei : cFaces)
        {
            label nbrCelli = -1;

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                // The nbrCell is actually the master (let it handle the face)
                if (nbrCelli <= newCelli)
                {
                    nbrCelli = -1;
                }
            }

            nbr.push_back(nbrCelli);
        }

        Foam::sortedOrder(nbr, order);

        for (const label index : order)
        {
            if (nbr[index] >= 0)
            {
                oldToNewFace[cFaces[index]] = newFacei++;
            }
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position" << " for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// Determine face order such that inside region faces are sorted
// upper-triangular but inbetween region faces are handled like boundary faces.
labelList getRegionFaceOrder
(
    const primitiveMesh& mesh,
    const labelList& cellOrder,     // New to old cell
    const labelList& cellToRegion   // Old cell to region
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));

    labelList oldToNewFace(mesh.nFaces(), -1);

    label newFacei = 0;

    DynamicList<label> nbr(64);
    DynamicList<label> order(64);

    forAll(cellOrder, newCelli)
    {
        const label oldCelli = cellOrder[newCelli];

        const cell& cFaces = mesh.cells()[oldCelli];

        // Neighbouring cells
        nbr.clear();

        for (const label facei : cFaces)
        {
            label nbrCelli = -1;

            if (mesh.isInternalFace(facei))
            {
                // Internal face. Get cell on other side.
                nbrCelli = reverseCellOrder[mesh.faceNeighbour()[facei]];
                if (nbrCelli == newCelli)
                {
                    nbrCelli = reverseCellOrder[mesh.faceOwner()[facei]];
                }

                // The nbrCell is actually the master (let it handle the face)
                if (nbrCelli <= newCelli)
                {
                    nbrCelli = -1;
                }
                else
                {
                    // A region boundary? - treat like an external face
                    label ownRegion = cellToRegion[oldCelli];
                    label neiRegion = cellToRegion[cellOrder[nbrCelli]];

                    if (ownRegion != neiRegion)
                    {
                        nbrCelli = -1;
                    }
                }
            }

            nbr.push_back(nbrCelli);
        }

        Foam::sortedOrder(nbr, order);

        for (const label index : order)
        {
            if (nbr[index] >= 0)
            {
                oldToNewFace[cFaces[index]] = newFacei++;
            }
        }
    }

    // This seems to be broken...

    // Do region interfaces
    {
        const label nRegions = max(cellToRegion)+1;

        // Sort in increasing region
        SortableList<label> sortKey(mesh.nInternalFaces(), labelMax);

        for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
        {
            label ownRegion = cellToRegion[mesh.faceOwner()[facei]];
            label neiRegion = cellToRegion[mesh.faceNeighbour()[facei]];

            if (ownRegion != neiRegion)
            {
                sortKey[facei] =
                    min(ownRegion, neiRegion)*nRegions
                   +max(ownRegion, neiRegion);
            }
        }

        sortKey.sort();

        // Extract.
        forAll(sortKey, i)
        {
            label key = sortKey[i];

            if (key == labelMax)
            {
                break;
            }

            oldToNewFace[sortKey.indices()[i]] = newFacei++;
        }
    }

    // Leave patch faces intact.
    for (label facei = newFacei; facei < mesh.nFaces(); facei++)
    {
        oldToNewFace[facei] = facei;
    }


    // Check done all faces.
    forAll(oldToNewFace, facei)
    {
        if (oldToNewFace[facei] == -1)
        {
            FatalErrorInFunction
                << "Did not determine new position for face " << facei
                << abort(FatalError);
        }
    }

    return invert(mesh.nFaces(), oldToNewFace);
}


// cellOrder: old cell for every new cell
// faceOrder: old face for every new face. Ordering of boundary faces not
//     changed.
autoPtr<mapPolyMesh> reorderMesh
(
    polyMesh& mesh,
    const labelList& cellOrder,
    const labelList& faceOrder
)
{
    labelList reverseCellOrder(invert(cellOrder.size(), cellOrder));
    labelList reverseFaceOrder(invert(faceOrder.size(), faceOrder));

    faceList newFaces(reorder(reverseFaceOrder, mesh.faces()));
    labelList newOwner
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceOwner())
        )
    );
    labelList newNeighbour
    (
        renumber
        (
            reverseCellOrder,
            reorder(reverseFaceOrder, mesh.faceNeighbour())
        )
    );

    // Check if any faces need swapping.
    labelHashSet flipFaceFlux(newOwner.size());
    forAll(newNeighbour, facei)
    {
        if (newOwner[facei] > newNeighbour[facei])
        {
            std::swap(newOwner[facei], newNeighbour[facei]);
            newFaces[facei].flip();
            flipFaceFlux.insert(facei);
        }
    }

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    labelList patchSizes(patches.size());
    labelList patchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    labelListList patchPointMap(patches.size());

    forAll(patches, patchi)
    {
        patchSizes[patchi] = patches[patchi].size();
        patchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
        patchPointMap[patchi] = identity(patches[patchi].nPoints());
    }

    mesh.resetPrimitives
    (
        autoPtr<pointField>(),  // <- null: leaves points untouched
        autoPtr<faceList>::New(std::move(newFaces)),
        autoPtr<labelList>::New(std::move(newOwner)),
        autoPtr<labelList>::New(std::move(newNeighbour)),
        patchSizes,
        patchStarts,
        true
    );


    // Re-do the faceZones
    {
        faceZoneMesh& faceZones = mesh.faceZones();
        faceZones.clearAddressing();
        forAll(faceZones, zoneI)
        {
            faceZone& fZone = faceZones[zoneI];
            labelList newAddressing(fZone.size());
            boolList newFlipMap(fZone.size());
            forAll(fZone, i)
            {
                label oldFacei = fZone[i];
                newAddressing[i] = reverseFaceOrder[oldFacei];
                if (flipFaceFlux.found(newAddressing[i]))
                {
                    newFlipMap[i] = !fZone.flipMap()[i];
                }
                else
                {
                    newFlipMap[i] = fZone.flipMap()[i];
                }
            }
            labelList newToOld(sortedOrder(newAddressing));
            fZone.resetAddressing
            (
                labelUIndList(newAddressing, newToOld)(),
                boolUIndList(newFlipMap, newToOld)()
            );
        }
    }
    // Re-do the cellZones
    {
        cellZoneMesh& cellZones = mesh.cellZones();
        cellZones.clearAddressing();
        forAll(cellZones, zoneI)
        {
            cellZones[zoneI] = labelUIndList
            (
                reverseCellOrder,
                cellZones[zoneI]
            )();
            Foam::sort(cellZones[zoneI]);
        }
    }


    return autoPtr<mapPolyMesh>::New
    (
        mesh,                       // const polyMesh& mesh,
        mesh.nPoints(),             // nOldPoints,
        mesh.nFaces(),              // nOldFaces,
        mesh.nCells(),              // nOldCells,
        identity(mesh.nPoints()),   // pointMap,
        List<objectMap>(),          // pointsFromPoints,
        faceOrder,                  // faceMap,
        List<objectMap>(),          // facesFromPoints,
        List<objectMap>(),          // facesFromEdges,
        List<objectMap>(),          // facesFromFaces,
        cellOrder,                  // cellMap,
        List<objectMap>(),          // cellsFromPoints,
        List<objectMap>(),          // cellsFromEdges,
        List<objectMap>(),          // cellsFromFaces,
        List<objectMap>(),          // cellsFromCells,
        identity(mesh.nPoints()),   // reversePointMap,
        reverseFaceOrder,           // reverseFaceMap,
        reverseCellOrder,           // reverseCellMap,
        flipFaceFlux,               // flipFaceFlux,
        patchPointMap,              // patchPointMap,
        labelListList(),            // pointZoneMap,
        labelListList(),            // faceZonePointMap,
        labelListList(),            // faceZoneFaceMap,
        labelListList(),            // cellZoneMap,
        pointField(),               // preMotionPoints,
        patchStarts,                // oldPatchStarts,
        oldPatchNMeshPoints,        // oldPatchNMeshPoints
        autoPtr<scalarField>()      // oldCellVolumes
    );
}


// Return new to old cell numbering, region-wise
CompactListList<label> regionRenumber
(
    const renumberMethod& method,
    const fvMesh& mesh,
    const labelList& cellToRegion,
    label nRegions = -1   // Number of regions or auto-detect
)
{
    // Info<< "Determining cell order:" << endl;

    if (nRegions < 0)
    {
        nRegions = Foam::max(cellToRegion)+1;
    }

    // Initially the per-region cell selection
    CompactListList<label> regionCellOrder
    (
        invertOneToManyCompact(nRegions, cellToRegion)
    );

    if (method.no_topology())
    {
        // Special case when renumberMesh is only used for decomposition.
        // - can skip generating the connectivity
        // - nonetheless calculate the order in case it is non-identity

        timer.resetTimeIncrement();

        forAll(regionCellOrder, regioni)
        {
            // Note: cellMap is identical to regionToCells[regioni]
            // since it is already sorted

            labelList subCellOrder =
                method.renumber(regionCellOrder[regioni].size());

            // Per region reordering (inplace but with SubList)
            regionCellOrder[regioni] =
                labelUIndList(regionCellOrder[regioni], subCellOrder)();
        }

        timings[TimingType::RENUMBER] += timer.timeIncrement();
    }
    else if (method.needs_mesh())
    {
        timer.resetTimeIncrement();

        forAll(regionCellOrder, regioni)
        {
            // Info<< "    region " << regioni << " starts at "
            //     << regionCellOrder.localStart(regioni) << nl;

            // No parallel communication
            const bool oldParRun = UPstream::parRun(false);

            // Get local sub-mesh
            fvMeshSubset subsetter(mesh, regionCellOrder[regioni]);

            // Note: cellMap will be identical to regionToCells[regioni]
            // (assuming they are properly sorted!)
            const labelList& cellMap = subsetter.cellMap();

            labelList subCellOrder = method.renumber(subsetter.subMesh());

            UPstream::parRun(oldParRun);  // Restore parallel state

            // Per region reordering
            regionCellOrder[regioni] = labelUIndList(cellMap, subCellOrder);
        }

        timings[TimingType::RENUMBER] += timer.timeIncrement();
    }
    else
    {
        timer.resetTimeIncrement();

        // Create adjacency matrix of the full mesh and subset subsequently.
        // This is more efficient than creating adjacency matrices of
        // sub-meshes.

        // No parallel communication
        const bool oldParRun = UPstream::parRun(false);

        // The local connectivity of the full (non-subsetted) mesh
        CompactListList<label> meshCellCells;
        globalMeshData::calcCellCells(mesh, meshCellCells);
        UPstream::parRun(oldParRun);  // Restore parallel state

        timings[TimingType::CELL_CELLS] += timer.timeIncrement();

        // For the respective subMesh selections
        bitSet subsetCells(mesh.nCells());

        forAll(regionCellOrder, regioni)
        {
            // Info<< "    region " << regioni << " starts at "
            //     << regionCellOrder.localStart(regioni) << nl;

            subsetCells = false;
            subsetCells.set(regionCellOrder[regioni]);

            // Connectivity of local sub-mesh
            labelList cellMap;
            CompactListList<label> subCellCells =
                subsetAdjacency(subsetCells, meshCellCells, cellMap);

            timings[TimingType::CELL_CELLS] += timer.timeIncrement();

            // No parallel communication
            const bool oldParRun = UPstream::parRun(false);

            labelList subCellOrder = method.renumber(subCellCells);

            UPstream::parRun(oldParRun);  // Restore parallel state

            // Per region reordering
            regionCellOrder[regioni] = labelUIndList(cellMap, subCellOrder);

            timings[TimingType::RENUMBER] += timer.timeIncrement();
        }
    }
    // Info<< endl;

    return regionCellOrder;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noFunctionObjects();           // Never use function objects

    timeSelector::addOptions_singleTime();  // Single-time options

    argList::addNote
    (
        "Renumber mesh cells to reduce the bandwidth. Use the -lib option or"
        " dictionary 'libs' entry to load additional libraries"
    );

    argList::addDryRunOption
    (
        "Test without writing. "
        "Changes -write-maps to write VTK output."
    );
    argList::addVerboseOption();

    argList::addOption("dict", "file", "Alternative renumberMeshDict");

    argList::addBoolOption
    (
        "frontWidth",
        "Calculate the RMS of the front-width"
    );

    argList::addBoolOption
    (
        "decompose",
        "Aggregate initially with a decomposition method (serial only)"
    );

    argList::addBoolOption
    (
        "write-maps",
        "Write renumber mappings"
    );

    argList::addBoolOption
    (
        "no-fields",
        "Suppress renumbering of fields (eg, when they are only uniform)"
    );

    argList::addBoolOption
    (
        "list-renumber",
        "List available renumbering methods"
    );

    argList::addOption
    (
        "renumber-method",
        "name",
        "Specify renumber method (default: CuthillMcKee) without dictionary"
    );

    argList::addOption
    (
        "renumber-coeffs",
        "string-content",
        "Specify renumber coefficients (dictionary content) as string. "
        "eg, 'reverse true;'"
    );

    #include "addAllRegionOptions.H"
    #include "addOverwriteOption.H"

    // -------------------------

    #include "setRootCase.H"

    {
        bool listOptions = false;

        if (args.found("list-renumber"))
        {
            listOptions = true;
            Info<< nl
                << "Available renumber methods:" << nl
                << "    "
                << flatOutput(renumberMethod::supportedMethods()) << nl
                << nl;
        }

        if (listOptions)
        {
            return 0;
        }
    }

    const bool dryrun = args.dryRun();

    const bool readDict = args.found("dict");
    const bool doDecompose = args.found("decompose");
    const bool overwrite = args.found("overwrite");
    const bool doFields = !args.found("no-fields");
    const bool doFrontWidth = args.found("frontWidth") && !doDecompose;

    word renumberMethodName;
    args.readIfPresent("renumber-method", renumberMethodName);

    if (doDecompose && UPstream::parRun())
    {
        FatalErrorIn(args.executable())
            << "Cannot use -decompose option in parallel ... giving up" << nl
            << exit(FatalError);
    }

    #include "createTime.H"

    // Allow override of decomposeParDict location
    fileName decompDictFile;
    if
    (
        args.readIfPresent("decomposeParDict", decompDictFile)
     && !decompDictFile.empty() && !decompDictFile.isAbsolute()
    )
    {
        decompDictFile = runTime.globalPath()/decompDictFile;
    }

    // Get region names
    #include "getAllRegionOptions.H"

    // Set time from specified time options, or force start from Time=0
    timeSelector::setTimeIfPresent(runTime, args, true);

    // Capture current time information for non-overwrite
    const Tuple2<instant, label> startTime
    (
        instant(runTime.value(), runTime.timeName()),
        runTime.timeIndex()
    );

    // Start/reset all timings
    timer.resetTime();
    timings = Foam::zero{};

    #include "createNamedMeshes.H"

    timings[TimingType::READ_MESH] += timer.timeIncrement();


    for (fvMesh& mesh : meshes)
    {
        // Reset time in case of multiple meshes and not overwrite
        if (!overwrite)
        {
            runTime.setTime(startTime.first(), startTime.second());
        }

        const word oldInstance = mesh.pointsInstance();

        label band;
        scalar profile;
        scalar sumSqrIntersect;
        getBand
        (
            doFrontWidth,
            mesh.nCells(),
            mesh.faceOwner(),
            mesh.faceNeighbour(),
            band,
            profile,
            sumSqrIntersect
        );

        reduce(band, maxOp<label>());
        reduce(profile, sumOp<scalar>());

        Info<< "Mesh " << mesh.name()
            << " size: " << mesh.globalData().nTotalCells() << nl
            << "Before renumbering" << nl
            << "    band           : " << band << nl
            << "    profile        : " << profile << nl;

        if (doFrontWidth)
        {
            reduce(sumSqrIntersect, sumOp<scalar>());
            scalar rmsFrontwidth = Foam::sqrt
            (
                sumSqrIntersect/mesh.globalData().nTotalCells()
            );

            Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
        }

        Info<< endl;

        bool sortCoupledFaceCells = false;
        bool writeMaps = args.found("write-maps");
        bool orderPoints = false;
        bool useRegionFaceOrder = false;
        label blockSize = 0;

        // Construct renumberMethod
        dictionary renumberDict;
        autoPtr<renumberMethod> renumberPtr;

        if (readDict)
        {
            const word dictName("renumberMeshDict");
            #include "setSystemMeshDictionaryIO.H"

            Info<< "Renumber according to " << dictIO.name() << nl << endl;

            renumberDict = IOdictionary::readContents(dictIO);

            renumberPtr = renumberMethod::New(renumberDict);

            // This options are likely orthogonal to -decompose mode
            if (!doDecompose)
            {
                sortCoupledFaceCells =
                    renumberDict.getOrDefault("sortCoupledFaceCells", false);

                if (sortCoupledFaceCells)
                {
                    Info<< "Sorting cells on coupled boundaries to be last."
                        << nl << endl;
                }

                blockSize = renumberDict.getOrDefault<label>("blockSize", 0);

                if (blockSize > 0)
                {
                    Info<< "Ordering cells into regions of size " << blockSize
                        << " (using decomposition);"
                        << " ordering faces into region-internal"
                        << " and region-external."
                        << nl << endl;
                }

                if (blockSize > 0)
                {
                    useRegionFaceOrder =
                        renumberDict.getOrDefault("regionFaceOrder", false);
                }
            }

            orderPoints = renumberDict.getOrDefault("orderPoints", false);
            if (orderPoints)
            {
                Info<< "Ordering points into internal and boundary points."
                    << nl << endl;
            }

            if
            (
                renumberDict.readIfPresent("writeMaps", writeMaps)
             && writeMaps
            )
            {
                Info<< "Writing renumber maps (new to old) to polyMesh."
                    << nl << endl;
            }
        }
        else
        {
            if (args.found("renumber-coeffs"))
            {
                ITstream is = args.lookup("renumber-coeffs");

                is >> renumberDict;

                Info<< "Specified renumber coefficients:" << nl
                    << renumberDict << nl;
            }

            if (!renumberMethodName.empty())
            {
                // Specify/override the "method" within dictionary
                renumberDict.set("method", renumberMethodName);
                renumberPtr = renumberMethod::New(renumberDict);
            }
            else if (renumberDict.found("method"))
            {
                // Use the "method" type within dictionary
                renumberPtr = renumberMethod::New(renumberDict);
            }
        }

        // Default renumbering method
        if (!renumberPtr)
        {
            renumberPtr.reset(new CuthillMcKeeRenumber(renumberDict));
            Info<< "Using renumber-method: " << renumberPtr().type()
                << " [default]" << endl;
        }
        else
        {
            Info<< "Using renumber-method: " << renumberPtr().type()
                << endl;
        }


        // Read parallel reconstruct maps
        labelIOList cellProcAddressing
        (
            IOobject
            (
                "cellProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                (dryrun ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
                IOobject::AUTO_WRITE
            ),
            labelList()
        );

        labelIOList faceProcAddressing
        (
            IOobject
            (
                "faceProcAddressing",
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                (dryrun ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
                IOobject::AUTO_WRITE
            ),
            labelList()
        );
        labelIOList pointProcAddressing
        (
            IOobject
            (
                "pointProcAddressing",
                mesh.pointsInstance(),
                polyMesh::meshSubDir,
                mesh,
                (dryrun ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
                IOobject::AUTO_WRITE
            ),
            labelList()
        );
        labelIOList boundaryProcAddressing
        (
            IOobject
            (
                "boundaryProcAddressing",
                mesh.pointsInstance(),
                polyMesh::meshSubDir,
                mesh,
                (dryrun ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
                IOobject::AUTO_WRITE
            ),
            labelList()
        );


        // List of stored objects to clear from mesh (after reading)
        DynamicList<regIOobject*> storedObjects;

        if (!dryrun && doFields)
        {
            Info<< nl << "Reading fields" << nl;

            timer.resetTimeIncrement();

            IOobjectList objects(mesh, runTime.timeName());
            storedObjects.reserve(objects.size());

            const predicates::always nameMatcher;

            // Read GeometricFields

            #undef  doLocalCode
            #define doLocalCode(FieldType)                                    \
            readFields<FieldType>(mesh, objects, nameMatcher, storedObjects);

            // Read volume fields
            doLocalCode(volScalarField);
            doLocalCode(volVectorField);
            doLocalCode(volSphericalTensorField);
            doLocalCode(volSymmTensorField);
            doLocalCode(volTensorField);

            // Read internal fields
            doLocalCode(volScalarField::Internal);
            doLocalCode(volVectorField::Internal);
            doLocalCode(volSphericalTensorField::Internal);
            doLocalCode(volSymmTensorField::Internal);
            doLocalCode(volTensorField::Internal);

            // Read surface fields
            doLocalCode(surfaceScalarField);
            doLocalCode(surfaceVectorField);
            doLocalCode(surfaceSphericalTensorField);
            doLocalCode(surfaceSymmTensorField);
            doLocalCode(surfaceTensorField);


            // Read point fields
            const pointMesh& pMesh = pointMesh::New(mesh);
            #undef  doLocalCode
            #define doLocalCode(FieldType)                                    \
            readFields<FieldType>(pMesh, objects, nameMatcher, storedObjects);

            doLocalCode(pointScalarField);
            doLocalCode(pointVectorField);
            doLocalCode(pointSphericalTensorField);
            doLocalCode(pointSymmTensorField);
            doLocalCode(pointTensorField);

            #undef doLocalCode

            timings[TimingType::READ_FIELDS] += timer.timeIncrement();

            // Write loaded fields when mesh.write() is called
            for (auto* fldptr : storedObjects)
            {
                fldptr->writeOpt(IOobject::AUTO_WRITE);
            }
        }


        // Read sets
        PtrList<cellSet> cellSets;
        PtrList<faceSet> faceSets;
        PtrList<pointSet> pointSets;

        if (!dryrun)
        {
            // Read sets
            IOobjectList objects(mesh, mesh.facesInstance(), "polyMesh/sets");
            ReadFields(objects, cellSets);
            ReadFields(objects, faceSets);
            ReadFields(objects, pointSets);
        }


        Info<< endl;

        // From renumbering:
        // - from new cell/face back to original cell/face
        labelList cellOrder;
        labelList faceOrder;

        if (blockSize > 0 && !doDecompose)
        {
            timer.resetTimeIncrement();

            // Renumbering in two phases. Should be done in one so mapping of
            // fields is done correctly!

            label nBlocks = mesh.nCells()/blockSize;
            Info<< "nBlocks   = " << nBlocks << endl;

            // Read decompositionMethod dictionary
            dictionary decomposeDict(renumberDict.subDict("blockCoeffs"));
            decomposeDict.set("numberOfSubdomains", nBlocks);

            // No parallel communication
            const bool oldParRun = UPstream::parRun(false);

            autoPtr<decompositionMethod> decomposePtr =
                decompositionMethod::New(decomposeDict);

            labelList cellToRegion
            (
                decomposePtr().decompose
                (
                    mesh,
                    mesh.cellCentres()
                )
            );

            UPstream::parRun(oldParRun);  // Restore parallel state

            timings[TimingType::DECOMPOSE] += timer.timeIncrement();

            // For debugging: write out region
            createScalarField
            (
                mesh,
                "cellDist",
                cellToRegion
            )().write();

            Info<< nl << "Written decomposition as volScalarField to "
                << "cellDist for use in postprocessing."
                << nl << endl;


            CompactListList<label> regionCellOrder =
                regionRenumber
                (
                    renumberPtr(),
                    mesh,
                    cellToRegion,
                    decomposePtr().nDomains()
                );

            cellOrder = regionCellOrder.values();

            // Determine new to old face order with new cell numbering
            if (useRegionFaceOrder)
            {
                faceOrder = getRegionFaceOrder(mesh, cellOrder, cellToRegion);
            }
            else
            {
                faceOrder = getFaceOrder(mesh, cellOrder);
            }
        }
        else
        {
            if (doDecompose)
            {
                // Two-step renumbering.
                // 1. decompose into regions (like decomposePar)
                // 2. renumber each sub-region

                timer.resetTimeIncrement();

                // Read decompositionMethod dictionary
                IOdictionary decomposeDict
                (
                    IOobject::selectIO
                    (
                        IOobject
                        (
                            decompositionModel::canonicalName,
                            runTime.time().system(),
                            mesh.regionName(),  // region (if non-default)
                            runTime,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE,
                            IOobject::NO_REGISTER
                        ),
                        decompDictFile
                    )
                );

                // No parallel communication
                const bool oldParRun = UPstream::parRun(false);

                autoPtr<decompositionMethod> decomposePtr
                    = decompositionMethod::New(decomposeDict);

                labelList cellToRegion
                (
                    decomposePtr().decompose
                    (
                        mesh,
                        mesh.cellCentres()
                    )
                );

                timings[TimingType::DECOMPOSE] += timer.timeIncrement();

                UPstream::parRun(oldParRun);  // Restore parallel state

                CompactListList<label> regionCellOrder =
                    regionRenumber
                    (
                        renumberPtr(),
                        mesh,
                        cellToRegion,
                        decomposePtr().nDomains()
                    );

                cellOrder = regionCellOrder.values();

                // HACK - retain partitioning information for possible use
                // at a later stage
                mesh.data().set
                (
                    "requested.partition-offsets",
                    regionCellOrder.offsets()
                );
            }
            else
            {
                // Determines sorted back to original cell ordering

                const auto& method = renumberPtr();

                timer.resetTimeIncrement();

                if (method.no_topology())
                {
                    cellOrder = method.renumber(mesh.nCells());
                }
                else
                {
                    cellOrder = method.renumber(mesh);
                }

                timings[TimingType::RENUMBER] += timer.timeIncrement();
            }


            if (sortCoupledFaceCells)
            {
                // Change order so all coupled patch faceCells are at the end.
                const polyBoundaryMesh& pbm = mesh.boundaryMesh();

                // Collect all boundary cells on coupled patches
                label nBndCells = 0;
                forAll(pbm, patchi)
                {
                    if (pbm[patchi].coupled())
                    {
                        nBndCells += pbm[patchi].size();
                    }
                }

                labelList reverseCellOrder = invert(mesh.nCells(), cellOrder);

                labelList bndCells(nBndCells);
                labelList bndCellMap(nBndCells);
                nBndCells = 0;
                forAll(pbm, patchi)
                {
                    if (pbm[patchi].coupled())
                    {
                        const labelUList& faceCells = pbm[patchi].faceCells();
                        forAll(faceCells, i)
                        {
                            label celli = faceCells[i];

                            if (reverseCellOrder[celli] != -1)
                            {
                                bndCells[nBndCells] = celli;
                                bndCellMap[nBndCells++] =
                                    reverseCellOrder[celli];
                                reverseCellOrder[celli] = -1;
                            }
                        }
                    }
                }
                bndCells.resize(nBndCells);
                bndCellMap.resize(nBndCells);

                // Sort
                labelList order(Foam::sortedOrder(bndCellMap));

                // Redo newReverseCellOrder
                labelList newReverseCellOrder(mesh.nCells(), -1);

                label sortedI = mesh.nCells();
                forAllReverse(order, i)
                {
                    label origCelli = bndCells[order[i]];
                    newReverseCellOrder[origCelli] = --sortedI;
                }

                Info<< "Ordered all " << nBndCells
                    << " cells with a coupled face"
                    << " to the end of the cell list, starting at " << sortedI
                    << endl;

                // Compact
                sortedI = 0;
                forAll(cellOrder, newCelli)
                {
                    label origCelli = cellOrder[newCelli];
                    if (newReverseCellOrder[origCelli] == -1)
                    {
                        newReverseCellOrder[origCelli] = sortedI++;
                    }
                }

                // Update sorted back to original (unsorted) map
                cellOrder = invert(mesh.nCells(), newReverseCellOrder);
            }


            // Determine new to old face order with new cell numbering
            faceOrder = getFaceOrder(mesh, cellOrder);
        }


        if (!overwrite)
        {
            ++runTime;
        }


        // Change the mesh.
        autoPtr<mapPolyMesh> map = reorderMesh(mesh, cellOrder, faceOrder);


        if (orderPoints)
        {
            polyTopoChange meshMod(mesh);
            autoPtr<mapPolyMesh> pointOrderMap = meshMod.changeMesh
            (
                mesh,
                false,      // inflate
                true,       // syncParallel
                false,      // orderCells
                orderPoints // orderPoints
            );

            // Combine point reordering into map.
            const_cast<labelList&>(map().pointMap()) = labelUIndList
            (
                map().pointMap(),
                pointOrderMap().pointMap()
            )();

            inplaceRenumber
            (
                pointOrderMap().reversePointMap(),
                const_cast<labelList&>(map().reversePointMap())
            );
        }


        // Update fields
        mesh.updateMesh(map());

        // Update proc maps
        if (cellProcAddressing.headerOk())
        {
            if (returnReduceAnd(cellProcAddressing.size() == mesh.nCells()))
            {
                Info<< "Renumbering processor cell decomposition map "
                    << cellProcAddressing.name() << endl;

                cellProcAddressing = labelList
                (
                    labelUIndList(cellProcAddressing, map().cellMap())
                );
            }
            else
            {
                Info<< "Not writing inconsistent processor cell decomposition"
                    << " map " << cellProcAddressing.filePath() << endl;
                cellProcAddressing.writeOpt(IOobject::NO_WRITE);
            }
        }
        else
        {
            cellProcAddressing.writeOpt(IOobject::NO_WRITE);
        }

        if (faceProcAddressing.headerOk())
        {
            if (returnReduceAnd(faceProcAddressing.size() == mesh.nFaces()))
            {
                Info<< "Renumbering processor face decomposition map "
                    << faceProcAddressing.name() << endl;

                faceProcAddressing = labelList
                (
                    labelUIndList(faceProcAddressing, map().faceMap())
                );

                // Detect any flips.
                const labelHashSet& fff = map().flipFaceFlux();
                for (const label facei : fff)
                {
                    label masterFacei = faceProcAddressing[facei];

                    faceProcAddressing[facei] = -masterFacei;

                    if (masterFacei == 0)
                    {
                        FatalErrorInFunction
                            << " masterFacei:" << masterFacei
                            << exit(FatalError);
                    }
                }
            }
            else
            {
                Info<< "Not writing inconsistent processor face decomposition"
                    << " map " << faceProcAddressing.filePath() << endl;
                faceProcAddressing.writeOpt(IOobject::NO_WRITE);
            }
        }
        else
        {
            faceProcAddressing.writeOpt(IOobject::NO_WRITE);
        }

        if (pointProcAddressing.headerOk())
        {
            if (returnReduceAnd(pointProcAddressing.size() == mesh.nPoints()))
            {
                Info<< "Renumbering processor point decomposition map "
                    << pointProcAddressing.name() << endl;

                pointProcAddressing = labelList
                (
                    labelUIndList(pointProcAddressing, map().pointMap())
                );
            }
            else
            {
                Info<< "Not writing inconsistent processor point decomposition"
                    << " map " << pointProcAddressing.filePath() << endl;
                pointProcAddressing.writeOpt(IOobject::NO_WRITE);
            }
        }
        else
        {
            pointProcAddressing.writeOpt(IOobject::NO_WRITE);
        }

        if (boundaryProcAddressing.headerOk())
        {
            if
            (
                returnReduceAnd
                (
                    boundaryProcAddressing.size() == mesh.boundaryMesh().size()
                )
            )
            {
                // No renumbering needed
            }
            else
            {
                Info<< "Not writing inconsistent processor patch decomposition"
                    << " map " << boundaryProcAddressing.filePath() << endl;
                boundaryProcAddressing.writeOpt(IOobject::NO_WRITE);
            }
        }
        else
        {
            boundaryProcAddressing.writeOpt(IOobject::NO_WRITE);
        }




        // Move mesh (since morphing might not do this)
        if (map().hasMotionPoints())
        {
            mesh.movePoints(map().preMotionPoints());
        }


        {
            label band;
            scalar profile;
            scalar sumSqrIntersect;
            getBand
            (
                doFrontWidth,
                mesh.nCells(),
                mesh.faceOwner(),
                mesh.faceNeighbour(),
                band,
                profile,
                sumSqrIntersect
            );
            reduce(band, maxOp<label>());
            reduce(profile, sumOp<scalar>());

            Info<< "After renumbering";
            if (doDecompose)
            {
                Info<< " [values are misleading with -decompose option]";
            }

            Info<< nl
                << "    band           : " << band << nl
                << "    profile        : " << profile << nl;

            if (doFrontWidth)
            {
                reduce(sumSqrIntersect, sumOp<scalar>());
                scalar rmsFrontwidth = Foam::sqrt
                (
                    sumSqrIntersect/mesh.globalData().nTotalCells()
                );

                Info<< "    rms frontwidth : " << rmsFrontwidth << nl;
            }

            Info<< endl;
        }

        if (orderPoints)
        {
            // Force edge calculation (since only reason that points would
            // need to be sorted)
            (void)mesh.edges();

            label nTotPoints = returnReduce
            (
                mesh.nPoints(),
                sumOp<label>()
            );
            label nTotIntPoints = returnReduce
            (
                mesh.nInternalPoints(),
                sumOp<label>()
            );

            label nTotEdges = returnReduce
            (
                mesh.nEdges(),
                sumOp<label>()
            );
            label nTotIntEdges = returnReduce
            (
                mesh.nInternalEdges(),
                sumOp<label>()
            );
            label nTotInt0Edges = returnReduce
            (
                mesh.nInternal0Edges(),
                sumOp<label>()
            );
            label nTotInt1Edges = returnReduce
            (
                mesh.nInternal1Edges(),
                sumOp<label>()
            );

            Info<< "Points:" << nl
                << "    total   : " << nTotPoints << nl
                << "    internal: " << nTotIntPoints << nl
                << "    boundary: " << nTotPoints-nTotIntPoints << nl
                << "Edges:" << nl
                << "    total   : " << nTotEdges << nl
                << "    internal: " << nTotIntEdges << nl
                << "        internal using 0 boundary points: "
                << nTotInt0Edges << nl
                << "        internal using 1 boundary points: "
                << nTotInt1Edges-nTotInt0Edges << nl
                << "        internal using 2 boundary points: "
                << nTotIntEdges-nTotInt1Edges << nl
                << "    boundary: " << nTotEdges-nTotIntEdges << nl
                << endl;
        }


        if (dryrun)
        {
            if (writeMaps)
            {
                fileName file = mesh.time().globalPath()/"renumberMesh";

                const word& regionDir = mesh.regionName();

                if (!regionDir.empty())
                {
                    file += "_" + regionDir;
                }

                // VTK output
                {
                    const vtk::vtuCells cells(mesh);

                    vtk::internalMeshWriter writer
                    (
                        mesh,
                        cells,
                        file,
                        UPstream::parRun()
                    );

                    writer.writeGeometry();
                    writer.beginCellData();
                    writer.writeCellIDs();

                    labelList ids;

                    if (UPstream::parRun())
                    {
                        const label cellOffset =
                            mesh.globalData().globalMeshCellAddr().localStart();

                        ids.resize(mesh.nCells());
                        std::transform
                        (
                            map().cellMap().cbegin(),
                            map().cellMap().cend(),
                            ids.begin(),
                            [=](const label id) { return (id + cellOffset); }
                        );

                        writer.writeCellData("origCellID", ids);
                        writer.writeProcIDs();
                    }
                    else
                    {
                        writer.writeCellData("origCellID", map().cellMap());

                        // Recover any partition information (in serial)
                        globalIndex partitionOffsets;
                        if
                        (
                            mesh.data().readIfPresent
                            (
                                "requested.partition-offsets",
                                partitionOffsets.offsets()
                            )
                        )
                        {
                            if (partitionOffsets.totalSize() != mesh.nCells())
                            {
                                WarningInFunction
                                    << "Requested partition total-size: "
                                    << partitionOffsets.totalSize()
                                    << " mesh total-size: "
                                    << mesh.nCells()
                                    << " ... ignoring" << endl;

                                partitionOffsets.clear();
                            }
                        }

                        ids.resize(partitionOffsets.totalSize());

                        for (const label proci : partitionOffsets.allProcs())
                        {
                            ids.slice(partitionOffsets.range(proci)) = proci;
                        }

                        if (!partitionOffsets.empty())
                        {
                            writer.writeCellData("procID", ids);
                        }
                    }

                    Info<< "Wrote renumbered mesh to "
                        << mesh.time().relativePath(writer.output())
                        << " for visualization."
                        << endl;
                }
            }
        }
        else
        {
            timer.resetTimeIncrement();

            if (overwrite)
            {
                mesh.setInstance(oldInstance);
            }
            else
            {
                mesh.setInstance(runTime.timeName());
            }

            Info<< "Writing mesh to " << mesh.facesInstance() << endl;

            // Remove old procAddressing files
            processorMeshes::removeFiles(mesh);

            // Update refinement data

            hexRef8Data refData
            (
                IOobject
                (
                    "dummy",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    (dryrun ? IOobject::NO_READ : IOobject::READ_IF_PRESENT),
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                )
            );
            refData.updateMesh(map());
            refData.write();

            // Update sets
            topoSet::updateMesh(mesh.facesInstance(), map(), cellSets);
            topoSet::updateMesh(mesh.facesInstance(), map(), faceSets);
            topoSet::updateMesh(mesh.facesInstance(), map(), pointSets);

            mesh.write();

            timings[TimingType::WRITING] += timer.timeIncrement();

            if (writeMaps)
            {
                // For debugging: write out region
                createScalarField
                (
                    mesh,
                    "origCellID",
                    map().cellMap()
                )().write();

                createScalarField
                (
                    mesh,
                    "cellID",
                    identity(mesh.nCells())
                )().write();

                Info<< nl
                    << "Wrote current cellID and origCellID as volScalarField"
                    << " for use in postprocessing." << nl << endl;

                IOobject meshMapIO
                (
                    "map-name",
                    mesh.facesInstance(),
                    polyMesh::meshSubDir,
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                );

                meshMapIO.resetHeader("cellMap");
                IOListRef<label>(meshMapIO, map().cellMap()).write();

                meshMapIO.resetHeader("faceMap");
                IOListRef<label>(meshMapIO, map().faceMap()).write();

                meshMapIO.resetHeader("pointMap");
                IOListRef<label>(meshMapIO, map().pointMap()).write();
            }
        }

        // Cleanup loaded fields
        while (!storedObjects.empty())
        {
            storedObjects.back()->checkOut();
            storedObjects.pop_back();
        }
    }

    Info<< nl
        << "Timings:" << nl
        << "    read mesh   : " << timings[TimingType::READ_MESH] << nl
        << "    read fields : " << timings[TimingType::READ_FIELDS] << nl
        << "    decompose   : " << timings[TimingType::DECOMPOSE] << nl
        << "    cell-cells  : " << timings[TimingType::CELL_CELLS] << nl
        << "    renumber    : " << timings[TimingType::RENUMBER] << nl
        << "    write       : " << timings[TimingType::WRITING] << nl
        << "TotalTime = " << timer.elapsedTime() << " s" << nl
        << nl;

    runTime.printExecutionTime(Info);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
