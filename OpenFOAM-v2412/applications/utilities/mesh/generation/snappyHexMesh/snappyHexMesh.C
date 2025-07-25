/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Application
    snappyHexMesh

Group
    grpMeshGenerationUtilities

Description
    Automatic split hex mesher. Refines and snaps to surface.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "fvMesh.H"
#include "snappyRefineDriver.H"
#include "snappySnapDriver.H"
#include "snappyLayerDriver.H"
#include "searchableSurfaces.H"
#include "refinementSurfaces.H"
#include "refinementFeatures.H"
#include "shellSurfaces.H"
#include "decompositionMethod.H"
#include "fvMeshDistribute.H"
#include "wallPolyPatch.H"
#include "refinementParameters.H"
#include "snapParameters.H"
#include "layerParameters.H"
#include "vtkCoordSetWriter.H"
#include "vtkSurfaceWriter.H"
#include "faceSet.H"
#include "motionSmoother.H"
#include "polyTopoChange.H"
#include "indirectPrimitivePatch.H"
#include "surfZoneIdentifierList.H"
#include "UnsortedMeshedSurface.H"
#include "MeshedSurface.H"
#include "globalIndex.H"
#include "IOmanip.H"
#include "decompositionModel.H"
#include "fvMeshTools.H"
#include "profiling.H"
#include "processorMeshes.H"
#include "snappyVoxelMeshDriver.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convert size (as fraction of defaultCellSize) to refinement level
label sizeCoeffToRefinement
(
    const scalar level0Coeff,   // ratio of hex cell size v.s. defaultCellSize
    const scalar sizeCoeff
)
{
     return round(::log(level0Coeff/sizeCoeff)/::log(2));
}


autoPtr<refinementSurfaces> createRefinementSurfaces
(
    const searchableSurfaces& allGeometry,
    const dictionary& surfacesDict,
    const dictionary& shapeControlDict,
    const label gapLevelIncrement,
    const scalar level0Coeff
)
{
    autoPtr<refinementSurfaces> surfacePtr;

    // Count number of surfaces.
    label surfi = 0;
    forAll(allGeometry.names(), geomi)
    {
        const word& geomName = allGeometry.names()[geomi];

        if (surfacesDict.found(geomName))
        {
            surfi++;
        }
    }

    labelList surfaces(surfi);
    wordList names(surfi);
    PtrList<surfaceZonesInfo> surfZones(surfi);

    labelList regionOffset(surfi);

    labelList globalMinLevel(surfi, Zero);
    labelList globalMaxLevel(surfi, Zero);
    labelList globalLevelIncr(surfi, Zero);
    PtrList<dictionary> globalPatchInfo(surfi);
    List<Map<label>> regionMinLevel(surfi);
    List<Map<label>> regionMaxLevel(surfi);
    List<Map<label>> regionLevelIncr(surfi);
    List<Map<scalar>> regionAngle(surfi);
    List<Map<autoPtr<dictionary>>> regionPatchInfo(surfi);

    wordHashSet unmatchedKeys(surfacesDict.toc());

    surfi = 0;
    forAll(allGeometry.names(), geomi)
    {
        const word& geomName = allGeometry.names()[geomi];

        const entry* ePtr = surfacesDict.findEntry(geomName, keyType::REGEX);

        if (ePtr)
        {
            const dictionary& shapeDict = ePtr->dict();
            unmatchedKeys.erase(ePtr->keyword());

            names[surfi] = geomName;
            surfaces[surfi] = geomi;

            const searchableSurface& surface = allGeometry[geomi];

            // Find the index in shapeControlDict
            // Invert surfaceCellSize to get the refinementLevel

            const word scsFuncName =
                shapeDict.get<word>("surfaceCellSizeFunction");

            const dictionary& scsDict =
                shapeDict.optionalSubDict(scsFuncName + "Coeffs");

            const scalar surfaceCellSize =
                scsDict.get<scalar>("surfaceCellSizeCoeff");

            const label refLevel = sizeCoeffToRefinement
            (
                level0Coeff,
                surfaceCellSize
            );

            globalMinLevel[surfi] = refLevel;
            globalMaxLevel[surfi] = refLevel;
            globalLevelIncr[surfi] = gapLevelIncrement;

            // Surface zones
            surfZones.set
            (
                surfi,
                new surfaceZonesInfo
                (
                    surface,
                    shapeDict,
                    allGeometry.regionNames()[surfaces[surfi]]
                )
            );


            // Global perpendicular angle
            if (shapeDict.found("patchInfo"))
            {
                globalPatchInfo.set
                (
                    surfi,
                    shapeDict.subDict("patchInfo").clone()
                );
            }


            // Per region override of patchInfo

            if (shapeDict.found("regions"))
            {
                const dictionary& regionsDict = shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfi]].regions();

                forAll(regionNames, regioni)
                {
                    if (regionsDict.found(regionNames[regioni]))
                    {
                        // Get the dictionary for region
                        const dictionary& regionDict = regionsDict.subDict
                        (
                            regionNames[regioni]
                        );

                        if (regionDict.found("patchInfo"))
                        {
                            regionPatchInfo[surfi].insert
                            (
                                regioni,
                                regionDict.subDict("patchInfo").clone()
                            );
                        }
                    }
                }
            }

            // Per region override of cellSize
            if (shapeDict.found("regions"))
            {
                const dictionary& shapeControlRegionsDict =
                    shapeDict.subDict("regions");
                const wordList& regionNames =
                    allGeometry[surfaces[surfi]].regions();

                forAll(regionNames, regioni)
                {
                    if (shapeControlRegionsDict.found(regionNames[regioni]))
                    {
                        const dictionary& shapeControlRegionDict =
                            shapeControlRegionsDict.subDict
                            (
                                regionNames[regioni]
                            );

                        const word scsFuncName =
                            shapeControlRegionDict.get<word>
                            (
                                "surfaceCellSizeFunction"
                            );
                        const dictionary& scsDict =
                            shapeControlRegionDict.subDict
                            (
                                scsFuncName + "Coeffs"
                            );

                        const scalar surfaceCellSize =
                            scsDict.get<scalar>("surfaceCellSizeCoeff");

                        const label refLevel = sizeCoeffToRefinement
                        (
                            level0Coeff,
                            surfaceCellSize
                        );

                        regionMinLevel[surfi].insert(regioni, refLevel);
                        regionMaxLevel[surfi].insert(regioni, refLevel);
                        regionLevelIncr[surfi].insert(regioni, 0);
                    }
                }
            }

            surfi++;
        }
    }

    // Calculate local to global region offset
    label nRegions = 0;

    forAll(surfaces, surfi)
    {
        regionOffset[surfi] = nRegions;
        nRegions += allGeometry[surfaces[surfi]].regions().size();
    }

    // Rework surface specific information into information per global region
    labelList minLevel(nRegions, Zero);
    labelList maxLevel(nRegions, Zero);
    labelList gapLevel(nRegions, -1);
    PtrList<dictionary> patchInfo(nRegions);

    forAll(globalMinLevel, surfi)
    {
        label nRegions = allGeometry[surfaces[surfi]].regions().size();

        // Initialise to global (i.e. per surface)
        for (label i = 0; i < nRegions; i++)
        {
            label globalRegioni = regionOffset[surfi] + i;
            minLevel[globalRegioni] = globalMinLevel[surfi];
            maxLevel[globalRegioni] = globalMaxLevel[surfi];
            gapLevel[globalRegioni] =
                maxLevel[globalRegioni]
              + globalLevelIncr[surfi];

            if (globalPatchInfo.set(surfi))
            {
                patchInfo.set
                (
                    globalRegioni,
                    globalPatchInfo[surfi].clone()
                );
            }
        }

        // Overwrite with region specific information
        forAllConstIters(regionMinLevel[surfi], iter)
        {
            label globalRegioni = regionOffset[surfi] + iter.key();

            minLevel[globalRegioni] = iter();
            maxLevel[globalRegioni] = regionMaxLevel[surfi][iter.key()];
            gapLevel[globalRegioni] =
                maxLevel[globalRegioni]
              + regionLevelIncr[surfi][iter.key()];
        }

        const Map<autoPtr<dictionary>>& localInfo = regionPatchInfo[surfi];
        forAllConstIters(localInfo, iter)
        {
            label globalRegioni = regionOffset[surfi] + iter.key();
            patchInfo.set(globalRegioni, iter()().clone());
        }
    }

    surfacePtr.reset
    (
        new refinementSurfaces
        (
            allGeometry,
            surfaces,
            names,
            surfZones,
            regionOffset,
            minLevel,
            maxLevel,
            gapLevel,
            scalarField(nRegions, -GREAT),  //perpendicularAngle,
            patchInfo,
            false                           //dryRun
        )
    );


    const refinementSurfaces& rf = surfacePtr();

    // Determine maximum region name length
    label maxLen = 0;
    forAll(rf.surfaces(), surfi)
    {
        label geomi = rf.surfaces()[surfi];
        const wordList& regionNames = allGeometry.regionNames()[geomi];
        forAll(regionNames, regioni)
        {
            maxLen = Foam::max(maxLen, label(regionNames[regioni].size()));
        }
    }


    Info<< setw(maxLen) << "Region"
        << setw(10) << "Min Level"
        << setw(10) << "Max Level"
        << setw(10) << "Gap Level" << nl
        << setw(maxLen) << "------"
        << setw(10) << "---------"
        << setw(10) << "---------"
        << setw(10) << "---------" << endl;

    forAll(rf.surfaces(), surfi)
    {
        label geomi = rf.surfaces()[surfi];

        Info<< rf.names()[surfi] << ':' << nl;

        const wordList& regionNames = allGeometry.regionNames()[geomi];

        forAll(regionNames, regioni)
        {
            label globali = rf.globalRegion(surfi, regioni);

            Info<< setw(maxLen) << regionNames[regioni]
                << setw(10) << rf.minLevel()[globali]
                << setw(10) << rf.maxLevel()[globali]
                << setw(10) << rf.gapLevel()[globali] << endl;
        }
    }


    return surfacePtr;
}


void extractSurface
(
    const polyMesh& mesh,
    const Time& runTime,
    const labelHashSet& includePatches,
    const fileName& outFileName
)
{
    const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

    // Collect sizes. Hash on names to handle local-only patches (e.g.
    //  processor patches)
    HashTable<label> patchSize(1024);
    label nFaces = 0;
    for (const label patchi : includePatches)
    {
        const polyPatch& pp = bMesh[patchi];
        patchSize.insert(pp.name(), pp.size());
        nFaces += pp.size();
    }
    Pstream::mapCombineGather(patchSize, plusEqOp<label>());


    // Allocate zone/patch for all patches
    HashTable<label> compactZoneID(1024);
    if (Pstream::master())
    {
        forAllConstIters(patchSize, iter)
        {
            compactZoneID.insert(iter.key(), compactZoneID.size());
        }
    }
    Pstream::broadcast(compactZoneID);


    // Rework HashTable into labelList just for speed of conversion
    labelList patchToCompactZone(bMesh.size(), -1);
    forAllConstIters(compactZoneID, iter)
    {
        label patchi = bMesh.findPatchID(iter.key());
        if (patchi != -1)
        {
            patchToCompactZone[patchi] = iter.val();
        }
    }

    // Collect faces on zones
    DynamicList<label> faceLabels(nFaces);
    DynamicList<label> compactZones(nFaces);
    for (const label patchi : includePatches)
    {
        const polyPatch& pp = bMesh[patchi];
        forAll(pp, i)
        {
            faceLabels.append(pp.start()+i);
            compactZones.append(patchToCompactZone[pp.index()]);
        }
    }

    // Addressing engine for all faces
    uindirectPrimitivePatch allBoundary
    (
        UIndirectList<face>(mesh.faces(), faceLabels),
        mesh.points()
    );


    // Find correspondence to master points
    labelList pointToGlobal;
    labelList uniqueMeshPoints;
    autoPtr<globalIndex> globalNumbers = mesh.globalData().mergePoints
    (
        allBoundary.meshPoints(),
        allBoundary.meshPointMap(),
        pointToGlobal,
        uniqueMeshPoints
    );

    // Gather all unique points on master
    List<pointField> gatheredPoints(Pstream::nProcs());
    gatheredPoints[Pstream::myProcNo()] = pointField
    (
        mesh.points(),
        uniqueMeshPoints
    );
    Pstream::gatherList(gatheredPoints);

    // Gather all faces
    List<faceList> gatheredFaces(Pstream::nProcs());
    gatheredFaces[Pstream::myProcNo()] = allBoundary.localFaces();
    forAll(gatheredFaces[Pstream::myProcNo()], i)
    {
        inplaceRenumber(pointToGlobal, gatheredFaces[Pstream::myProcNo()][i]);
    }
    Pstream::gatherList(gatheredFaces);

    // Gather all ZoneIDs
    List<labelList> gatheredZones(Pstream::nProcs());
    gatheredZones[Pstream::myProcNo()].transfer(compactZones);
    Pstream::gatherList(gatheredZones);

    // On master combine all points, faces, zones
    if (Pstream::master())
    {
        pointField allPoints = ListListOps::combine<pointField>
        (
            gatheredPoints,
            accessOp<pointField>()
        );
        gatheredPoints.clear();

        faceList allFaces = ListListOps::combine<faceList>
        (
            gatheredFaces,
            accessOp<faceList>()
        );
        gatheredFaces.clear();

        labelList allZones = ListListOps::combine<labelList>
        (
            gatheredZones,
            accessOp<labelList>()
        );
        gatheredZones.clear();


        // Zones
        surfZoneIdentifierList surfZones(compactZoneID.size());
        forAllConstIters(compactZoneID, iter)
        {
            surfZones[iter()] = surfZoneIdentifier(iter.key(), iter());
            Info<< "surfZone " << iter()  <<  " : " << surfZones[iter()].name()
                << endl;
        }

        UnsortedMeshedSurface<face> unsortedFace
        (
            std::move(allPoints),
            std::move(allFaces),
            std::move(allZones),
            surfZones
        );


        MeshedSurface<face> sortedFace(unsortedFace);

        fileName globalCasePath
        (
            runTime.processorCase()
          ? runTime.globalPath()/outFileName
          : runTime.path()/outFileName
        );
        globalCasePath.clean();  // Remove unneeded ".."

        Info<< "Writing merged surface to " << globalCasePath << endl;

        sortedFace.write(globalCasePath);
    }
}


label checkAlignment(const polyMesh& mesh, const scalar tol, Ostream& os)
{
    // Check all edges aligned with one of the coordinate axes
    const faceList& faces = mesh.faces();
    const pointField& points = mesh.points();

    label nUnaligned = 0;

    forAll(faces, facei)
    {
        const face& f = faces[facei];
        forAll(f, fp)
        {
            label fp1 = f.fcIndex(fp);
            const linePointRef e(edge(f[fp], f[fp1]).line(points));
            const vector v(e.vec());
            const scalar magV(mag(v));
            if (magV > ROOTVSMALL)
            {
                for
                (
                    direction dir = 0;
                    dir < pTraits<vector>::nComponents;
                    ++dir
                )
                {
                    const scalar s(mag(v[dir]));
                    if (s > magV*tol && s < magV*(1-tol))
                    {
                        ++nUnaligned;
                        break;
                    }
                }
            }
        }
    }

    reduce(nUnaligned, sumOp<label>());

    if (nUnaligned)
    {
        os << "Initial mesh has " << nUnaligned
            << " edges unaligned with any of the coordinate axes" << nl << endl;
    }
    return nUnaligned;
}


// Check writing tolerance before doing any serious work
scalar getMergeDistance
(
    const polyMesh& mesh,
    const scalar mergeTol,
    const bool dryRun
)
{
    const boundBox& meshBb = mesh.bounds();
    scalar mergeDist = mergeTol * meshBb.mag();

    Info<< nl
        << "Overall mesh bounding box  : " << meshBb << nl
        << "Relative tolerance         : " << mergeTol << nl
        << "Absolute matching distance : " << mergeDist << nl
        << endl;

    // check writing tolerance
    if (mesh.time().writeFormat() == IOstreamOption::ASCII && !dryRun)
    {
        const scalar writeTol = std::pow
        (
            scalar(10),
            -scalar(IOstream::defaultPrecision())
        );

        if (mergeTol < writeTol)
        {
            FatalErrorInFunction
                << "Your current settings specify ASCII writing with "
                << IOstream::defaultPrecision() << " digits precision." << nl
                << "Your merging tolerance (" << mergeTol
                << ") is finer than this." << nl
                << "Change to binary writeFormat, "
                << "or increase the writePrecision" << endl
                << "or adjust the merge tolerance (mergeTol)."
                << exit(FatalError);
        }
    }

    return mergeDist;
}


void removeZeroSizedPatches(fvMesh& mesh)
{
    // Remove any zero-sized ones. Assumes
    // - processor patches are already only there if needed
    // - all other patches are available on all processors
    // - but coupled ones might still be needed, even if zero-size
    //   (e.g. processorCyclic)
    // See also logic in createPatch.
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    labelList oldToNew(pbm.size(), -1);
    label newPatchi = 0;
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if
            (
                isA<coupledPolyPatch>(pp)
             || returnReduceOr(pp.size())
            )
            {
                // Coupled (and unknown size) or uncoupled and used
                oldToNew[patchi] = newPatchi++;
            }
        }
    }

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            oldToNew[patchi] = newPatchi++;
        }
    }


    const label nKeepPatches = newPatchi;

    // Shuffle unused ones to end
    if (nKeepPatches != pbm.size())
    {
        Info<< endl
            << "Removing zero-sized patches:" << endl << incrIndent;

        forAll(oldToNew, patchi)
        {
            if (oldToNew[patchi] == -1)
            {
                Info<< indent << pbm[patchi].name()
                    << " type " << pbm[patchi].type()
                    << " at position " << patchi << endl;
                oldToNew[patchi] = newPatchi++;
            }
        }
        Info<< decrIndent;

        fvMeshTools::reorderPatches(mesh, oldToNew, nKeepPatches, true);
        Info<< endl;
    }
}


// Write mesh and additional information
void writeMesh
(
    const string& msg,
    const meshRefinement& meshRefiner,
    const meshRefinement::debugType debugLevel,
    const meshRefinement::writeType writeLevel
)
{
    // Note: don't want to use time().cpuTimeIncrement since layer addition
    // somehow resets timer ...
    cpuTime timer;

    const fvMesh& mesh = meshRefiner.mesh();

    meshRefiner.printMeshInfo(debugLevel, msg, true);
    Info<< "Writing mesh to time " << meshRefiner.timeName() << endl;

    processorMeshes::removeFiles(mesh);
    if (!debugLevel && !(writeLevel&meshRefinement::WRITELAYERSETS))
    {
        topoSet::removeFiles(mesh);
    }
    refinementHistory::removeFiles(mesh);

    meshRefiner.write
    (
        debugLevel,
        meshRefinement::writeType(writeLevel | meshRefinement::WRITEMESH),
        mesh.time().path()/meshRefiner.timeName()
    );
    Info<< "Wrote mesh in = "
        << timer.cpuTimeIncrement() << " s." << endl;
}


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Automatic split hex mesher. Refines and snaps to surface"
    );

    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    #include "addProfilingOption.H"
    argList::addBoolOption
    (
        "checkGeometry",
        "Check all surface geometry for quality"
    );
    argList::addDryRunOption
    (
        "Check case set-up only using a single time step"
    );
    argList::addOption
    (
        "surfaceSimplify",
        "boundBox",
        "Simplify the surface using snappyHexMesh starting from a boundBox"
    );
    argList::addOption
    (
        "patches",
        "(patch0 .. patchN)",
        "Only triangulate selected patches (wildcards supported)"
    );
    argList::addOption
    (
        "outFile",
        "file",
        "Name of the file to save the simplified surface to"
    );
    argList::addOption("dict", "file", "Alternative snappyHexMeshDict");

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"

    const bool overwrite = args.found("overwrite");
    const bool checkGeometry = args.found("checkGeometry");
    const bool surfaceSimplify = args.found("surfaceSimplify");
    const bool dryRun = args.dryRun();

    if (dryRun)
    {
        Info<< "Operating in dry-run mode to detect set-up errors"
            << nl << endl;
    }


    #include "createNamedMesh.H"
    Info<< "Read mesh in = "
        << runTime.cpuTimeIncrement() << " s" << endl;

    // Check patches and faceZones are synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    meshRefinement::checkCoupledFaceZones(mesh);

    if (dryRun)
    {
        // Check if mesh aligned with cartesian axes
        checkAlignment(mesh, 1e-6, Pout);   //FatalIOError);
    }



    // Read meshing dictionary
    const word dictName("snappyHexMeshDict");
    #include "setSystemMeshDictionaryIO.H"
    const IOdictionary meshDict(dictIO);

    // Overall mesh generation mode
    const meshRefinement::MeshType meshType
    (
        meshRefinement::MeshTypeNames.getOrDefault
        (
            "type",
            meshDict,
            meshRefinement::CASTELLATED
        )
    );


    // all surface geometry
    const dictionary& geometryDict =
        meshRefinement::subDict(meshDict, "geometry", dryRun);

    // refinement parameters
    const dictionary& refineDict =
        meshRefinement::subDict(meshDict, "castellatedMeshControls", dryRun);

    // mesh motion and mesh quality parameters
    const dictionary& motionDict =
        meshRefinement::subDict(meshDict, "meshQualityControls", dryRun);

    // snap-to-surface parameters
    const dictionary& snapDict =
        meshRefinement::subDict(meshDict, "snapControls", dryRun);

    // layer addition parameters
    const dictionary& layerDict =
        meshRefinement::subDict(meshDict, "addLayersControls", dryRun);

    // absolute merge distance
    const scalar mergeDist = getMergeDistance
    (
        mesh,
        meshRefinement::get<scalar>
        (
            meshDict,
            "mergeTolerance",
            dryRun
        ),
        dryRun
    );

    const bool keepPatches(meshDict.getOrDefault("keepPatches", false));

    // Writer for writing lines
    autoPtr<coordSetWriter> setFormatter;
    {
        const word writerType
        (
            meshDict.getOrDefault<word>
            (
                "setFormat",
                coordSetWriters::vtkWriter::typeName    // Default: "vtk"
            )
        );

        setFormatter = coordSetWriter::New
        (
            writerType,
            meshDict.subOrEmptyDict("formatOptions").optionalSubDict(writerType)
        );
    }
    // Writer for writing surfaces
    refPtr<surfaceWriter> surfFormatter;
    {
        const word type
        (
            meshDict.getOrDefault<word>
            (
                "surfaceFormat",
                surfaceWriters::vtkWriter::typeName      // Default: "vtk"
            )
        );
        surfFormatter = surfaceWriter::New
        (
            type,
            meshDict.subOrEmptyDict("formatOptions").optionalSubDict(type)
        );
    }


    const scalar maxSizeRatio
    (
        meshDict.getOrDefault<scalar>("maxSizeRatio", 100)
    );


    // Read decomposePar dictionary
    dictionary decomposeDict;
    if (Pstream::parRun())
    {
        // Ensure demand-driven decompositionMethod finds alternative
        // decomposeParDict location properly.

        IOdictionary* dictPtr = new IOdictionary
        (
            IOobject::selectIO
            (
                IOobject
                (
                    decompositionModel::canonicalName,
                    runTime.system(),
                    runTime,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    IOobject::REGISTER
                ),
                args.getOrDefault<fileName>("decomposeParDict", "")
            )
        );

        // Store it on the object registry, but to be found it must also
        // have the expected "decomposeParDict" name.

        dictPtr->rename(decompositionModel::canonicalName);
        runTime.store(dictPtr);

        decomposeDict = *dictPtr;
    }
    else
    {
        decomposeDict.add("method", "none");
        decomposeDict.add("numberOfSubdomains", 1);
    }


    // Debug
    // ~~~~~

    // Set debug level
    meshRefinement::debugType debugLevel = meshRefinement::debugType
    (
        meshDict.getOrDefault<label>
        (
            "debug",
            0
        )
    );
    {
        wordList flags;
        if (meshDict.readIfPresent("debugFlags", flags))
        {
            debugLevel = meshRefinement::debugType
            (
                meshRefinement::readFlags
                (
                    meshRefinement::debugTypeNames,
                    flags
                )
            );
        }
    }
    if (debugLevel > 0)
    {
        meshRefinement::debug   = debugLevel;
        snappyRefineDriver::debug = debugLevel;
        snappySnapDriver::debug   = debugLevel;
        snappyLayerDriver::debug  = debugLevel;
    }

    // Set file writing level
    {
        wordList flags;
        if (meshDict.readIfPresent("writeFlags", flags))
        {
            meshRefinement::writeLevel
            (
                meshRefinement::writeType
                (
                    meshRefinement::readFlags
                    (
                        meshRefinement::writeTypeNames,
                        flags
                    )
                )
            );
        }
    }

    //// Set output level
    //{
    //    wordList flags;
    //    if (meshDict.readIfPresent("outputFlags", flags))
    //    {
    //        meshRefinement::outputLevel
    //        (
    //            meshRefinement::outputType
    //            (
    //                meshRefinement::readFlags
    //                (
    //                    meshRefinement::outputTypeNames,
    //                    flags
    //                )
    //            )
    //        );
    //    }
    //}

    // for the impatient who want to see some output files:
    profiling::writeNow();

    // Read geometry
    // ~~~~~~~~~~~~~

    searchableSurfaces allGeometry
    (
        IOobject
        (
            "abc",                      // dummy name
            mesh.time().constant(),     // instance
            //mesh.time().findInstance("triSurface", word::null),// instance
            "triSurface",               // local
            mesh.time(),                // registry
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        geometryDict,
        meshDict.getOrDefault("singleRegionName", true)
    );


    // Read refinement surfaces
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    autoPtr<refinementSurfaces> surfacesPtr;

    Info<< "Reading refinement surfaces." << endl;

    if (surfaceSimplify)
    {
        addProfiling(surfaceSimplify, "snappyHexMesh::surfaceSimplify");
        IOdictionary foamyHexMeshDict
        (
           IOobject
           (
                "foamyHexMeshDict",
                runTime.system(),
                runTime,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
           )
        );

        const dictionary& conformationDict =
            foamyHexMeshDict.subDict("surfaceConformation").subDict
            (
                "geometryToConformTo"
            );

        const dictionary& motionDict =
            foamyHexMeshDict.subDict("motionControl");

        const dictionary& shapeControlDict =
            motionDict.subDict("shapeControlFunctions");

        // Calculate current ratio of hex cells v.s. wanted cell size
        const scalar defaultCellSize =
            motionDict.get<scalar>("defaultCellSize");

        const scalar initialCellSize = ::pow(mesh.V()[0], 1.0/3.0);

        //Info<< "Wanted cell size  = " << defaultCellSize << endl;
        //Info<< "Current cell size = " << initialCellSize << endl;
        //Info<< "Fraction          = " << initialCellSize/defaultCellSize
        //    << endl;

        surfacesPtr =
            createRefinementSurfaces
            (
                allGeometry,
                conformationDict,
                shapeControlDict,
                refineDict.getOrDefault("gapLevelIncrement", 0),
                initialCellSize/defaultCellSize
            );

        profiling::writeNow();
    }
    else
    {
        surfacesPtr.reset
        (
            new refinementSurfaces
            (
                allGeometry,
                meshRefinement::subDict
                (
                    refineDict,
                    "refinementSurfaces",
                    dryRun
                ),
                refineDict.getOrDefault("gapLevelIncrement", 0),
                dryRun
            )
        );

        Info<< "Read refinement surfaces in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    refinementSurfaces& surfaces = surfacesPtr();


    // Checking only?

    if (checkGeometry)
    {
        // Check geometry amongst itself (e.g. intersection, size differences)

        // Extract patchInfo
        List<wordList> patchTypes(allGeometry.size());

        const PtrList<dictionary>& patchInfo = surfaces.patchInfo();
        const labelList& surfaceGeometry = surfaces.surfaces();
        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];
            const wordList& regNames = allGeometry.regionNames()[geomi];

            patchTypes[geomi].setSize(regNames.size());
            forAll(regNames, regioni)
            {
                label globalRegioni = surfaces.globalRegion(surfi, regioni);

                if (patchInfo.set(globalRegioni))
                {
                    patchTypes[geomi][regioni] =
                        meshRefinement::get<word>
                        (
                            patchInfo[globalRegioni],
                            "type",
                            dryRun,
                            keyType::REGEX,
                            word::null
                        );
                }
                else
                {
                    patchTypes[geomi][regioni] = wallPolyPatch::typeName;
                }
            }
        }

        // Write some stats
        allGeometry.writeStats(patchTypes, Info);
        // Check topology
        allGeometry.checkTopology(true);
        // Check geometry
        allGeometry.checkGeometry
        (
            maxSizeRatio,   // max size ratio
            1e-9,           // intersection tolerance
            setFormatter,
            0.01,           // min triangle quality
            true
        );

        if (!dryRun)
        {
            return 0;
        }
    }


    if (dryRun)
    {
        // Check geometry to mesh bounding box
        Info<< "Checking for geometry size relative to mesh." << endl;
        const boundBox& meshBb = mesh.bounds();
        forAll(allGeometry, geomi)
        {
            const searchableSurface& s = allGeometry[geomi];
            const boundBox& bb = s.bounds();

            scalar ratio = bb.mag() / meshBb.mag();
            if (ratio > maxSizeRatio || ratio < 1.0/maxSizeRatio)
            {
                Warning
                    << "    " << allGeometry.names()[geomi]
                    << " bounds differ from mesh"
                    << " by more than a factor " << maxSizeRatio << ":" << nl
                    << "        bounding box      : " << bb << nl
                    << "        mesh bounding box : " << meshBb
                    << endl;
            }
            if (!meshBb.contains(bb))
            {
                Warning
                    << "    " << allGeometry.names()[geomi]
                    << " bounds not fully contained in mesh" << nl
                    << "        bounding box      : " << bb << nl
                    << "        mesh bounding box : " << meshBb
                    << endl;
            }
        }
        Info<< endl;
    }




    // Read refinement shells
    // ~~~~~~~~~~~~~~~~~~~~~~

    Info<< "Reading refinement shells." << endl;
    shellSurfaces shells
    (
        allGeometry,
        meshRefinement::subDict(refineDict, "refinementRegions", dryRun),
        dryRun
    );
    Info<< "Read refinement shells in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    Info<< "Setting refinement level of surface to be consistent"
        << " with shells." << endl;
    surfaces.setMinLevelFields(shells);
    Info<< "Checked shell refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    // Optionally read limit shells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const dictionary limitDict(refineDict.subOrEmptyDict("limitRegions"));

    if (!limitDict.empty())
    {
        Info<< "Reading limit shells." << endl;
    }

    shellSurfaces limitShells(allGeometry, limitDict, dryRun);

    if (!limitDict.empty())
    {
        Info<< "Read limit shells in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    if (dryRun)
    {
        // Check for use of all geometry
        const wordList& allGeomNames = allGeometry.names();

        labelHashSet unusedGeometries(identity(allGeomNames.size()));
        unusedGeometries.erase(surfaces.surfaces());
        unusedGeometries.erase(shells.shells());
        unusedGeometries.erase(limitShells.shells());

        if (unusedGeometries.size())
        {
            IOWarningInFunction(geometryDict)
                << "The following geometry entries are not used:" << nl;
            for (const label geomi : unusedGeometries)
            {
                Info<< "    " << allGeomNames[geomi] << nl;
            }
            Info<< endl;
        }
    }




    // Read feature meshes
    // ~~~~~~~~~~~~~~~~~~~

    Info<< "Reading features." << endl;
    refinementFeatures features
    (
        mesh,
        PtrList<dictionary>
        (
            meshRefinement::lookup(refineDict, "features", dryRun)
        ),
        dryRun
    );
    Info<< "Read features in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;


    if (dryRun)
    {
        // Check geometry to mesh bounding box
        Info<< "Checking for line geometry size relative to surface geometry."
            << endl;

        OStringStream os;
        bool hasErrors = features.checkSizes
        (
            maxSizeRatio,   //const scalar maxRatio,
            mesh.bounds(),
            true,           //const bool report,
            os              //FatalIOError
        );
        if (hasErrors)
        {
            Warning<< os.str() << endl;
        }
    }


    // Refinement engine
    // ~~~~~~~~~~~~~~~~~

    Info<< nl
        << "Determining initial surface intersections" << nl
        << "-----------------------------------------" << nl
        << endl;

    // Main refinement engine
    meshRefinement meshRefiner
    (
        mesh,
        mergeDist,          // tolerance used in sorting coordinates
        overwrite,          // overwrite mesh files?
        surfaces,           // for surface intersection refinement
        features,           // for feature edges/point based refinement
        shells,             // for volume (inside/outside) refinement
        limitShells,        // limit of volume refinement
        labelList(),        // initial faces to test
        meshType,           // how to operate
        dryRun
    );

    if (!dryRun)
    {
        meshRefiner.updateIntersections(identity(mesh.nFaces()));
        Info<< "Calculated surface intersections in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }

    // Some stats
    meshRefiner.printMeshInfo(debugLevel, "Initial mesh", true);

    meshRefiner.write
    (
        meshRefinement::debugType(debugLevel&meshRefinement::OBJINTERSECTIONS),
        meshRefinement::writeType(0),
        mesh.time().path()/meshRefiner.timeName()
    );


    // Refinement parameters
    const refinementParameters refineParams(refineDict, dryRun);

    // Snap parameters
    const snapParameters snapParams(snapDict, dryRun);


    Info<< "Setting refinement level of surface to be consistent"
        << " with curvature." << endl;
    surfaces.setCurvatureMinLevelFields
    (
        refineParams.curvature(),
        meshRefiner.meshCutter().level0EdgeLength()
    );
    Info<< "Checked curvature refinement in = "
        << mesh.time().cpuTimeIncrement() << " s" << nl << endl;



    // Add all the cellZones and faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // 1. cellZones relating to surface (faceZones added later)

    const labelList namedSurfaces
    (
        surfaceZonesInfo::getNamedSurfaces(surfaces.surfZones())
    );

    labelList surfaceToCellZone = surfaceZonesInfo::addCellZonesToMesh
    (
        surfaces.surfZones(),
        namedSurfaces,
        mesh
    );


    // 2. cellZones relating to locations

    refineParams.addCellZonesToMesh(mesh);



    // Add all the surface regions as patches
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    //- Global surface region to patch (non faceZone surface) or patches
    //  (faceZone surfaces)
    labelList globalToMasterPatch;
    labelList globalToSlavePatch;


    {
        Info<< nl
            << "Adding patches for surface regions" << nl
            << "----------------------------------" << nl
            << endl;

        // From global region number to mesh patch.
        globalToMasterPatch.setSize(surfaces.nRegions(), -1);
        globalToSlavePatch.setSize(surfaces.nRegions(), -1);

        if (!dryRun)
        {
            Info<< setf(ios_base::left)
                << setw(6) << "Patch"
                << setw(20) << "Type"
                << setw(30) << "Region" << nl
                << setw(6) << "-----"
                << setw(20) << "----"
                << setw(30) << "------" << endl;
        }

        const labelList& surfaceGeometry = surfaces.surfaces();
        const PtrList<dictionary>& surfacePatchInfo = surfaces.patchInfo();
        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        forAll(surfaceGeometry, surfi)
        {
            label geomi = surfaceGeometry[surfi];

            const wordList& regNames = allGeometry.regionNames()[geomi];

            if (!dryRun)
            {
                Info<< surfaces.names()[surfi] << ':' << nl << nl;
            }

            const wordList& fzNames =
                surfaces.surfZones()[surfi].faceZoneNames();

            if (fzNames.empty())
            {
                // 'Normal' surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    label patchi;

                    if (surfacePatchInfo.set(globalRegioni))
                    {
                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            surfacePatchInfo[globalRegioni]
                        );
                    }
                    else
                    {
                        dictionary patchInfo;
                        patchInfo.set("type", wallPolyPatch::typeName);

                        patchi = meshRefiner.addMeshedPatch
                        (
                            regNames[i],
                            patchInfo
                        );
                    }

                    if (!dryRun)
                    {
                        Info<< setf(ios_base::left)
                            << setw(6) << patchi
                            << setw(20) << pbm[patchi].type()
                            << setw(30) << regNames[i] << nl;
                    }

                    globalToMasterPatch[globalRegioni] = patchi;
                    globalToSlavePatch[globalRegioni] = patchi;
                }
            }
            else
            {
                // Zoned surface
                forAll(regNames, i)
                {
                    label globalRegioni = surfaces.globalRegion(surfi, i);

                    // Add master side patch
                    {
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                regNames[i],
                                patchInfo
                            );
                        }

                        if (!dryRun)
                        {
                            Info<< setf(ios_base::left)
                                << setw(6) << patchi
                                << setw(20) << pbm[patchi].type()
                                << setw(30) << regNames[i] << nl;
                        }

                        globalToMasterPatch[globalRegioni] = patchi;
                    }
                    // Add slave side patch
                    {
                        const word slaveName = regNames[i] + "_slave";
                        label patchi;

                        if (surfacePatchInfo.set(globalRegioni))
                        {
                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                surfacePatchInfo[globalRegioni]
                            );
                        }
                        else
                        {
                            dictionary patchInfo;
                            patchInfo.set("type", wallPolyPatch::typeName);

                            patchi = meshRefiner.addMeshedPatch
                            (
                                slaveName,
                                patchInfo
                            );
                        }

                        if (!dryRun)
                        {
                            Info<< setf(ios_base::left)
                                << setw(6) << patchi
                                << setw(20) << pbm[patchi].type()
                                << setw(30) << slaveName << nl;
                        }

                        globalToSlavePatch[globalRegioni] = patchi;
                    }
                }

                // For now: have single faceZone per surface. Use first
                // region in surface for patch for zoning
                if (regNames.size())
                {
                    forAll(fzNames, fzi)
                    {
                        const word& fzName = fzNames[fzi];
                        label globalRegioni = surfaces.globalRegion(surfi, fzi);

                        meshRefiner.addFaceZone
                        (
                            fzName,
                            pbm[globalToMasterPatch[globalRegioni]].name(),
                            pbm[globalToSlavePatch[globalRegioni]].name(),
                            surfaces.surfZones()[surfi].faceType()
                        );
                    }
                }
            }

            if (!dryRun)
            {
                Info<< nl;
            }
        }
        Info<< "Added patches in = "
            << mesh.time().cpuTimeIncrement() << " s" << nl << endl;
    }



    // Add all information for all the remaining faceZones
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    HashTable<Pair<word>> faceZoneToPatches;
    forAll(mesh.faceZones(), zonei)
    {
        const word& fzName = mesh.faceZones()[zonei].name();

        label mpI, spI;
        surfaceZonesInfo::faceZoneType fzType;
        bool hasInfo = meshRefiner.getFaceZoneInfo(fzName, mpI, spI, fzType);

        if (!hasInfo)
        {
            // faceZone does not originate from a surface but presumably
            // from a cellZone pair instead
            string::size_type i = fzName.find("_to_");
            if (i != string::npos)
            {
                word cz0 = fzName.substr(0, i);
                word cz1 = fzName.substr(i+4, fzName.size()-i+4);
                word slaveName(cz1 + "_to_" + cz0);
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
            else
            {
                // Add as fzName + fzName_slave
                const word slaveName = fzName + "_slave";
                faceZoneToPatches.insert(fzName, Pair<word>(fzName, slaveName));
            }
        }
    }

    if (faceZoneToPatches.size())
    {
        snappyRefineDriver::addFaceZones
        (
            meshRefiner,
            refineParams,
            faceZoneToPatches
        );
    }



    // Re-do intersections on meshed boundaries since they use an extrapolated
    // other side
    {
        const labelList adaptPatchIDs(meshRefiner.meshedPatches());

        const polyBoundaryMesh& pbm = mesh.boundaryMesh();

        label nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            nFaces += pbm[adaptPatchIDs[i]].size();
        }

        labelList faceLabels(nFaces);
        nFaces = 0;
        forAll(adaptPatchIDs, i)
        {
            const polyPatch& pp = pbm[adaptPatchIDs[i]];
            forAll(pp, i)
            {
                faceLabels[nFaces++] = pp.start()+i;
            }
        }
        meshRefiner.updateIntersections(faceLabels);
    }



    // Parallel
    // ~~~~~~~~

    // Construct decomposition engine. Note: cannot use decompositionModel
    // MeshObject since we're clearing out the mesh inside the mesh generation.
    autoPtr<decompositionMethod> decomposerPtr
    (
        decompositionMethod::New
        (
            decomposeDict
        )
    );
    decompositionMethod& decomposer = *decomposerPtr;

    if (Pstream::parRun() && !decomposer.parallelAware())
    {
        FatalErrorInFunction
            << "You have selected decomposition method "
            << decomposer.typeName
            << " which is not parallel aware." << endl
            << "Please select one that is (hierarchical, ptscotch)"
            << exit(FatalError);
    }

    // Mesh distribution engine (uses tolerance to reconstruct meshes)
    fvMeshDistribute distributor(mesh);





    // Now do the real work -refinement -snapping -layers
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    const bool wantRefine
    (
        meshRefinement::get<bool>(meshDict, "castellatedMesh", dryRun)
    );
    const bool wantSnap
    (
        meshRefinement::get<bool>(meshDict, "snap", dryRun)
    );
    const bool wantLayers
    (
        meshRefinement::get<bool>(meshDict, "addLayers", dryRun)
    );

    if (dryRun)
    {
        string errorMsg(FatalError.message());
        string IOerrorMsg(FatalIOError.message());

        if (errorMsg.size() || IOerrorMsg.size())
        {
            //errorMsg = "[dryRun] " + errorMsg;
            //errorMsg.replaceAll("\n", "\n[dryRun] ");
            //IOerrorMsg = "[dryRun] " + IOerrorMsg;
            //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

            Warning
                << nl
                << "Missing/incorrect required dictionary entries:" << nl
                << nl
                << IOerrorMsg.c_str() << nl
                << errorMsg.c_str() << nl << nl
                << "Exiting dry-run" << nl << endl;

            FatalError.clear();
            FatalIOError.clear();

            return 0;
        }
    }


    // How to treat co-planar faces
    meshRefinement::FaceMergeType mergeType =
        meshRefinement::FaceMergeType::GEOMETRIC;
    {
        const bool mergePatchFaces
        (
            meshDict.getOrDefault("mergePatchFaces", true)
        );

        if (!mergePatchFaces)
        {
            Info<< "Not merging patch-faces of cell to preserve"
                << " (split)hex cell shape."
                << nl << endl;
            mergeType = meshRefinement::FaceMergeType::NONE;
        }
        else
        {
            const bool mergeAcrossPatches
            (
                meshDict.getOrDefault("mergeAcrossPatches", false)
            );

            if (mergeAcrossPatches)
            {
                Info<< "Merging co-planar patch-faces of cells"
                    << ", regardless of patch assignment"
                    << nl << endl;
                mergeType = meshRefinement::FaceMergeType::IGNOREPATCH;
            }
        }
    }



    if (wantRefine)
    {
        cpuTime timer;

        snappyRefineDriver refineDriver
        (
            meshRefiner,
            decomposer,
            distributor,
            globalToMasterPatch,
            globalToSlavePatch,
            setFormatter(),
            surfFormatter,
            dryRun
        );


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }


        refineDriver.doRefine
        (
            refineDict,
            refineParams,
            snapParams,
            refineParams.handleSnapProblems(),
            mergeType,
            motionDict
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches && !wantSnap && !wantLayers)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        if (!dryRun)
        {
            writeMesh
            (
                "Refined mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel()
            );
        }

        Info<< "Mesh refined in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }

    if (wantSnap)
    {
        cpuTime timer;

        snappySnapDriver snapDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch,
            dryRun
        );

        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        // Use the resolveFeatureAngle from the refinement parameters
        scalar curvature = refineParams.curvature();
        scalar planarAngle = refineParams.planarAngle();

        snapDriver.doSnap
        (
            snapDict,
            motionDict,
            mergeType,
            curvature,
            planarAngle,
            snapParams
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches && !wantLayers)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        if (!dryRun)
        {
            writeMesh
            (
                "Snapped mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel()
            );
        }

        Info<< "Mesh snapped in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }

    if (wantLayers)
    {
        cpuTime timer;

        // Layer addition parameters
        const layerParameters layerParams
        (
            layerDict,
            mesh.boundaryMesh(),
            dryRun
        );

        snappyLayerDriver layerDriver
        (
            meshRefiner,
            globalToMasterPatch,
            globalToSlavePatch,
            dryRun
        );

        // Use the maxLocalCells from the refinement parameters
        const bool preBalance =
            returnReduceOr(mesh.nCells() >= refineParams.maxLocalCells());


        if (!overwrite && !debugLevel)
        {
            const_cast<Time&>(mesh.time())++;
        }

        layerDriver.doLayers
        (
            layerDict,
            motionDict,
            layerParams,
            mergeType,
            preBalance,
            decomposer,
            distributor
        );

        // Remove zero sized patches originating from faceZones
        if (!keepPatches)
        {
            fvMeshTools::removeEmptyPatches(mesh, true);
        }

        if (!dryRun)
        {
            writeMesh
            (
                "Layer mesh",
                meshRefiner,
                debugLevel,
                meshRefinement::writeLevel()
            );
        }

        Info<< "Layers added in = "
            << timer.cpuTimeIncrement() << " s." << endl;

        profiling::writeNow();
    }


    {
        addProfiling(checkMesh, "snappyHexMesh::checkMesh");

        // Check final mesh
        Info<< "Checking final mesh ..." << endl;
        faceSet wrongFaces(mesh, "wrongFaces", mesh.nFaces()/100);
        motionSmoother::checkMesh(false, mesh, motionDict, wrongFaces, dryRun);
        const label nErrors = returnReduce
        (
            wrongFaces.size(),
            sumOp<label>()
        );

        if (nErrors > 0)
        {
            Info<< "Finished meshing with " << nErrors << " illegal faces"
                << " (concave, zero area or negative cell pyramid volume)"
                << endl;
            wrongFaces.write();
        }
        else
        {
            Info<< "Finished meshing without any errors" << endl;
        }

        profiling::writeNow();
    }


    if (surfaceSimplify)
    {
        addProfiling(surfaceSimplify, "snappyHexMesh::surfaceSimplify");

        const polyBoundaryMesh& bMesh = mesh.boundaryMesh();

        labelHashSet includePatches(bMesh.size());

        if (args.found("patches"))
        {
            includePatches = bMesh.patchSet
            (
                args.getList<wordRe>("patches")
            );
        }
        else
        {
            forAll(bMesh, patchi)
            {
                const polyPatch& patch = bMesh[patchi];

                if (!isA<processorPolyPatch>(patch))
                {
                    includePatches.insert(patchi);
                }
            }
        }

        fileName outFileName
        (
            args.getOrDefault<fileName>
            (
                "outFile",
                "constant/triSurface/simplifiedSurface.stl"
            )
        );

        extractSurface
        (
            mesh,
            runTime,
            includePatches,
            outFileName
        );

        pointIOField cellCentres
        (
            IOobject
            (
                "internalCellCentres",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh.cellCentres()
        );

        cellCentres.write();
    }

    profiling::writeNow();

    Info<< "Finished meshing in = "
        << runTime.elapsedCpuTime() << " s." << endl;


    if (dryRun)
    {
        string errorMsg(FatalError.message());
        string IOerrorMsg(FatalIOError.message());

        if (errorMsg.size() || IOerrorMsg.size())
        {
            //errorMsg = "[dryRun] " + errorMsg;
            //errorMsg.replaceAll("\n", "\n[dryRun] ");
            //IOerrorMsg = "[dryRun] " + IOerrorMsg;
            //IOerrorMsg.replaceAll("\n", "\n[dryRun] ");

            Perr<< nl
                << "Missing/incorrect required dictionary entries:" << nl
                << nl
                << IOerrorMsg.c_str() << nl
                << errorMsg.c_str() << nl << nl
                << "Exiting dry-run" << nl << endl;

            FatalError.clear();
            FatalIOError.clear();

            return 0;
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
