/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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
    splitMeshRegions

Group
    grpMeshManipulationUtilities

Description
    Splits mesh into multiple regions.

    Each region is defined as a domain whose cells can all be reached by
    cell-face-cell walking without crossing
    - boundary faces
    - additional faces from faceset (-blockedFaces faceSet).
    - any face between differing cellZones (-cellZones)

    Output is:
    - volScalarField with regions as different scalars (-detectOnly)
            or
    - mesh with multiple regions and mapped patches. These patches
      either cover the whole interface between two region (default) or
      only part according to faceZones (-useFaceZones)
            or
    - mesh with cells put into cellZones (-makeCellZones)

    Note:

    - multiple cellZones can be combined into a single region (cluster)
    for further analysis using the 'addZones' or 'combineZones' option:
        -addZones '((allSolids zoneA "zoneB.*")(allFluids none otherZone))'
    or
        -combineZones '((zoneA "zoneB.*")(none otherZone))
    This can be combined with e.g. 'cellZones' or 'cellZonesOnly'. The
    addZones option supplies the destination region name as first element in
    the list. The combineZones option synthesises the region name e.g.
        zoneA_zoneB0_zoneB1

    - cellZonesOnly does not do a walk and uses the cellZones only. Use
    this if you don't mind having disconnected domains in a single region.
    This option requires all cells to be in one (and one only) cellZone.

    - cellZonesFileOnly behaves like -cellZonesOnly but reads the cellZones
    from the specified file. This allows one to explicitly specify the region
    distribution and still have multiple cellZones per region.

    - prefixRegion prefixes all normal patches with region name (interface
    (patches already have region name prefix)

    - Should work in parallel.
    cellZones can differ on either side of processor boundaries in which case
    the faces get moved from processor patch to mapped patch. Not very well
    tested.

    - If a cell zone gets split into more than one region it can detect
    the largest matching region (-sloppyCellZones). This will accept any
    region that covers more than 50% of the zone. It has to be a subset
    so cannot have any cells in any other zone.

    - If explicitly a single region has been selected (-largestOnly or
      -insidePoint) its region name will be either
        - name of a cellZone it matches to or
        - "largestOnly" respectively "insidePoint" or
        - polyMesh::defaultRegion if additionally -overwrite
          (so it will overwrite the input mesh!)

    - writes maps like decomposePar back to original mesh:
        - pointRegionAddressing : for every point in this region the point in
        the original mesh
        - cellRegionAddressing  :   ,,      cell                ,,  cell    ,,
        - faceRegionAddressing  :   ,,      face                ,,  face in
        the original mesh + 'turning index'. For a face in the same orientation
        this is the original facelabel+1, for a turned face this is -facelabel-1
        - boundaryRegionAddressing : for every patch in this region the
        patch in the original mesh (or -1 if added patch)

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "regionSplit.H"
#include "fvMeshSubset.H"
#include "IOobjectList.H"
#include "volFields.H"
#include "faceSet.H"
#include "cellSet.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "edgeHashes.H"
#include "syncTools.H"
#include "ReadFields.H"
#include "mappedWallPolyPatch.H"
#include "fvMeshTools.H"
#include "processorMeshes.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Prepend prefix to selected patches.
void renamePatches
(
    fvMesh& mesh,
    const word& prefix,
    const labelList& patchesToRename
)
{
    polyBoundaryMesh& polyPatches =
        const_cast<polyBoundaryMesh&>(mesh.boundaryMesh());
    forAll(patchesToRename, i)
    {
        label patchi = patchesToRename[i];
        polyPatch& pp = polyPatches[patchi];

        if (isA<coupledPolyPatch>(pp))
        {
            WarningInFunction
                << "Encountered coupled patch " << pp.name()
                << ". Will only rename the patch itself,"
                << " not any referred patches."
                << " This might have to be done by hand."
                << endl;
        }

        pp.name() = prefix + '_' + pp.name();
    }
    // Recalculate any demand driven data (e.g. group to name lookup)
    polyPatches.updateMesh();
}


template<class GeoField>
void subsetVolFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& cellMap,
    const labelList& faceMap,
    const labelHashSet& addedPatches
)
{
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    for (const GeoField& fld : mesh.objectRegistry::csorted<GeoField>())
    {
        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] == Zero;
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt(IOobject::AUTO_WRITE);
        subFld->store();
    }
}


template<class GeoField>
void subsetSurfaceFields
(
    const fvMesh& mesh,
    const fvMesh& subMesh,
    const labelList& cellMap,
    const labelList& faceMap,
    const labelHashSet& addedPatches
)
{
    const labelList patchMap(identity(mesh.boundaryMesh().size()));

    for (const GeoField& fld : mesh.objectRegistry::csorted<GeoField>())
    {
        Info<< "Mapping field " << fld.name() << endl;

        tmp<GeoField> tSubFld
        (
            fvMeshSubset::interpolate
            (
                fld,
                subMesh,
                patchMap,
                cellMap,
                faceMap
            )
        );

        // Hack: set value to 0 for introduced patches (since don't
        //       get initialised.
        forAll(tSubFld().boundaryField(), patchi)
        {
            if (addedPatches.found(patchi))
            {
                tSubFld.ref().boundaryFieldRef()[patchi] == Zero;
            }
        }

        // Store on subMesh
        GeoField* subFld = tSubFld.ptr();
        subFld->rename(fld.name());
        subFld->writeOpt(IOobject::AUTO_WRITE);
        subFld->store();
    }
}

// Select all cells not in the region
labelList getNonRegionCells(const labelList& cellRegion, const label regionI)
{
    DynamicList<label> nonRegionCells(cellRegion.size());
    forAll(cellRegion, celli)
    {
        if (cellRegion[celli] != regionI)
        {
            nonRegionCells.append(celli);
        }
    }
    return nonRegionCells.shrink();
}


void addToInterface
(
    const polyMesh& mesh,
    const label zoneID,
    const label ownRegion,
    const label neiRegion,
    EdgeMap<Map<label>>& regionsToSize
)
{
    edge interface
    (
        min(ownRegion, neiRegion),
        max(ownRegion, neiRegion)
    );

    auto iter = regionsToSize.find(interface);

    if (iter.good())
    {
        // Check if zone present
        auto zoneIter = iter().find(zoneID);
        if (zoneIter.good())
        {
            // Found zone. Increment count.
            ++(*zoneIter);
        }
        else
        {
            // New or no zone. Insert with count 1.
            iter().insert(zoneID, 1);
        }
    }
    else
    {
        // Create new interface of size 1.
        regionsToSize(interface).insert(zoneID, 1);
    }
}


// Get region-region interface name and sizes.
// Returns interfaces as straight list for looping in identical order.
void getInterfaceSizes
(
    const polyMesh& mesh,
    const bool useFaceZones,
    const labelList& cellRegion,
    const wordList& regionNames,

    edgeList& interfaces,
    List<Pair<word>>& interfaceNames,
    labelList& interfaceSizes,
    labelList& faceToInterface
)
{
    // From region-region to faceZone (or -1) to number of faces.

    EdgeMap<Map<label>> regionsToSize;


    // Internal faces
    // ~~~~~~~~~~~~~~

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }

    // Boundary faces
    // ~~~~~~~~~~~~~~

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nBoundaryFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);

    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            addToInterface
            (
                mesh,
                (useFaceZones ? mesh.faceZones().whichZone(facei) : -1),
                ownRegion,
                neiRegion,
                regionsToSize
            );
        }
    }


    if (UPstream::parRun())
    {
        if (UPstream::master())
        {
            // Receive and add to my sizes
            for (const int proci : UPstream::subProcs())
            {
                EdgeMap<Map<label>> slaveSizes;
                IPstream::recv(slaveSizes, proci);

                forAllConstIters(slaveSizes, slaveIter)
                {
                    const Map<label>& slaveInfo = *slaveIter;

                    auto masterIter = regionsToSize.find(slaveIter.key());

                    if (masterIter.good())
                    {
                        // Same inter-region
                        Map<label>& masterInfo = *masterIter;

                        forAllConstIters(slaveInfo, iter)
                        {
                            const label zoneID = iter.key();
                            const label slaveSize = iter.val();

                            auto zoneIter = masterInfo.find(zoneID);
                            if (zoneIter.good())
                            {
                                *zoneIter += slaveSize;
                            }
                            else
                            {
                                masterInfo.insert(zoneID, slaveSize);
                            }
                        }
                    }
                    else
                    {
                        regionsToSize.insert(slaveIter.key(), slaveInfo);
                    }
                }
            }
        }
        else
        {
            // send to master
            OPstream::send(regionsToSize, UPstream::masterNo());
        }
    }

    // Rework

    Pstream::broadcast(regionsToSize);


    // Now we have the global sizes of all inter-regions.
    // Invert this on master and distribute.
    label nInterfaces = 0;
    forAllConstIters(regionsToSize, iter)
    {
        const Map<label>& info = iter.val();
        nInterfaces += info.size();
    }

    interfaces.setSize(nInterfaces);
    interfaceNames.setSize(nInterfaces);
    interfaceSizes.setSize(nInterfaces);
    EdgeMap<Map<label>> regionsToInterface(nInterfaces);

    nInterfaces = 0;
    forAllConstIters(regionsToSize, iter)
    {
        const edge& e = iter.key();
        const Map<label>& info = iter.val();

        const word& name0 = regionNames[e[0]];
        const word& name1 = regionNames[e[1]];

        forAllConstIters(info, infoIter)
        {
            interfaces[nInterfaces] = iter.key();
            label zoneID = infoIter.key();
            if (zoneID == -1)
            {
                interfaceNames[nInterfaces] = Pair<word>
                (
                    name0 + "_to_" + name1,
                    name1 + "_to_" + name0
                );
            }
            else
            {
                const word& zoneName = mesh.faceZones()[zoneID].name();
                interfaceNames[nInterfaces] = Pair<word>
                (
                    zoneName + "_" + name0 + "_to_" + name1,
                    zoneName + "_" + name1 + "_to_" + name0
                );
            }

            interfaceSizes[nInterfaces] = infoIter();
            regionsToInterface(e).insert(zoneID, nInterfaces);

            nInterfaces++;
        }
    }


    // Consistent interface information for all processors
    Pstream::broadcasts
    (
        UPstream::worldComm,
        interfaces,
        interfaceNames,
        interfaceSizes,
        regionsToInterface
    );

    // Mark all inter-region faces.
    faceToInterface.setSize(mesh.nFaces(), -1);

    forAll(mesh.faceNeighbour(), facei)
    {
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = cellRegion[mesh.faceNeighbour()[facei]];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
    forAll(coupledRegion, i)
    {
        label facei = i+mesh.nInternalFaces();
        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = coupledRegion[i];

        if (ownRegion != neiRegion)
        {
            label zoneID = -1;
            if (useFaceZones)
            {
                zoneID = mesh.faceZones().whichZone(facei);
            }

            edge interface
            (
                min(ownRegion, neiRegion),
                max(ownRegion, neiRegion)
            );

            faceToInterface[facei] = regionsToInterface[interface][zoneID];
        }
    }
}


// Create mesh for region.
autoPtr<mapPolyMesh> createRegionMesh
(
    const fvMesh& mesh,
    // Region info
    const labelList& cellRegion,
    const label regionI,
    const word& regionName,
    // Interface info
    const labelList& interfacePatches,
    const labelList& faceToInterface,

    autoPtr<fvMesh>& newMesh
)
{
    // Create dummy system/fv*
    fvMeshTools::createDummyFvMeshFiles(mesh, regionName, true);

    // Neighbour cellRegion.
    labelList coupledRegion(mesh.nBoundaryFaces());

    forAll(coupledRegion, i)
    {
        label celli = mesh.faceOwner()[i+mesh.nInternalFaces()];
        coupledRegion[i] = cellRegion[celli];
    }
    syncTools::swapBoundaryFaceList(mesh, coupledRegion);


    // Topology change container. Start off from existing mesh.
    polyTopoChange meshMod(mesh);

    // Cell remover engine
    removeCells cellRemover(mesh);

    // Select all but region cells
    labelList cellsToRemove(getNonRegionCells(cellRegion, regionI));

    // Find out which faces will get exposed. Note that this
    // gets faces in mesh face order. So both regions will get same
    // face in same order (important!)
    labelList exposedFaces = cellRemover.getExposedFaces(cellsToRemove);

    labelList exposedPatchIDs(exposedFaces.size());
    forAll(exposedFaces, i)
    {
        label facei = exposedFaces[i];
        label interfacei = faceToInterface[facei];

        label ownRegion = cellRegion[mesh.faceOwner()[facei]];
        label neiRegion = -1;

        if (mesh.isInternalFace(facei))
        {
            neiRegion = cellRegion[mesh.faceNeighbour()[facei]];
        }
        else
        {
            neiRegion = coupledRegion[facei-mesh.nInternalFaces()];
        }


        // Check which side is being kept - determines which of the two
        // patches will be used.

        label otherRegion = -1;

        if (ownRegion == regionI && neiRegion != regionI)
        {
            otherRegion = neiRegion;
        }
        else if (ownRegion != regionI && neiRegion == regionI)
        {
            otherRegion = ownRegion;
        }
        else
        {
            FatalErrorInFunction
                << "Exposed face:" << facei
                << " fc:" << mesh.faceCentres()[facei]
                << " has owner region " << ownRegion
                << " and neighbour region " << neiRegion
                << " when handling region:" << regionI
                << exit(FatalError);
        }

        // Find the patch.
        if (regionI < otherRegion)
        {
            exposedPatchIDs[i] = interfacePatches[interfacei];
        }
        else
        {
            exposedPatchIDs[i] = interfacePatches[interfacei]+1;
        }
    }

    // Remove faces
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        exposedPatchIDs,
        meshMod
    );

    autoPtr<mapPolyMesh> map = meshMod.makeMesh
    (
        newMesh,
        IOobject
        (
            regionName,
            mesh.time().timeName(),
            mesh.time(),
            IOobject::READ_IF_PRESENT,  // read fv* if present
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    return map;
}


void createAndWriteRegion
(
    const fvMesh& mesh,
    const labelList& cellRegion,
    const wordList& regionNames,
    const bool prefixRegion,
    const labelList& faceToInterface,
    const labelList& interfacePatches,
    const label regionI,
    const word& newMeshInstance
)
{
    Info<< "Creating mesh for region " << regionI
        << ' ' << regionNames[regionI] << endl;

    autoPtr<fvMesh> newMesh;
    autoPtr<mapPolyMesh> map = createRegionMesh
    (
        mesh,
        cellRegion,
        regionI,
        regionNames[regionI],
        interfacePatches,
        faceToInterface,
        newMesh
    );


    // Make map of all added patches
    labelHashSet addedPatches(2*interfacePatches.size());
    forAll(interfacePatches, interfacei)
    {
        addedPatches.insert(interfacePatches[interfacei]);
        addedPatches.insert(interfacePatches[interfacei]+1);
    }


    Info<< "Mapping fields" << endl;

    // Map existing fields
    newMesh().updateMesh(map());

    // Add subsetted fields
    subsetVolFields<volScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetVolFields<volTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );

    subsetSurfaceFields<surfaceScalarField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceVectorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSphericalTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceSymmTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );
    subsetSurfaceFields<surfaceTensorField>
    (
        mesh,
        newMesh(),
        map().cellMap(),
        map().faceMap(),
        addedPatches
    );


    const polyBoundaryMesh& newPatches = newMesh().boundaryMesh();
    newPatches.checkParallelSync(true);

    // Delete empty patches
    // ~~~~~~~~~~~~~~~~~~~~

    // Create reordering list to move patches-to-be-deleted to end
    labelList oldToNew(newPatches.size(), -1);
    DynamicList<label> sharedPatches(newPatches.size());
    label newI = 0;

    Info<< "Deleting empty patches" << endl;

    // Assumes all non-proc boundaries are on all processors!
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            if (returnReduceOr(pp.size()))
            {
                oldToNew[patchi] = newI;
                if (!addedPatches.found(patchi))
                {
                    sharedPatches.append(newI);
                }
                newI++;
            }
        }
    }

    // Same for processor patches (but need no reduction)
    forAll(newPatches, patchi)
    {
        const polyPatch& pp = newPatches[patchi];

        if (isA<processorPolyPatch>(pp) && pp.size())
        {
            oldToNew[patchi] = newI++;
        }
    }

    const label nNewPatches = newI;

    // Move all deleteable patches to the end
    forAll(oldToNew, patchi)
    {
        if (oldToNew[patchi] == -1)
        {
            oldToNew[patchi] = newI++;
        }
    }

    //reorderPatches(newMesh(), oldToNew, nNewPatches);
    fvMeshTools::reorderPatches(newMesh(), oldToNew, nNewPatches, true);

    // Rename shared patches with region name
    if (prefixRegion)
    {
        Info<< "Prefixing patches with region name" << endl;

        renamePatches(newMesh(), regionNames[regionI], sharedPatches);
    }


    Info<< "Writing new mesh" << endl;

    newMesh().setInstance(newMeshInstance);
    newMesh().write();
    topoSet::removeFiles(newMesh());
    processorMeshes::removeFiles(newMesh());

    // Write addressing files like decomposePar
    Info<< "Writing addressing to base mesh" << endl;

    labelIOList pointProcAddressing
    (
        IOobject
        (
            "pointRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        map().pointMap()
    );
    Info<< "Writing map " << pointProcAddressing.name()
        << " from region" << regionI
        << " points back to base mesh." << endl;
    pointProcAddressing.write();

    labelIOList faceProcAddressing
    (
        IOobject
        (
            "faceRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        newMesh().nFaces()
    );
    forAll(faceProcAddressing, facei)
    {
        // face + turning index. (see decomposePar)
        // Is the face pointing in the same direction?
        label oldFacei = map().faceMap()[facei];

        if
        (
            map().cellMap()[newMesh().faceOwner()[facei]]
         == mesh.faceOwner()[oldFacei]
        )
        {
            faceProcAddressing[facei] = oldFacei+1;
        }
        else
        {
            faceProcAddressing[facei] = -(oldFacei+1);
        }
    }
    Info<< "Writing map " << faceProcAddressing.name()
        << " from region" << regionI
        << " faces back to base mesh." << endl;
    faceProcAddressing.write();

    labelIOList cellProcAddressing
    (
        IOobject
        (
            "cellRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        map().cellMap()
    );
    Info<< "Writing map " <<cellProcAddressing.name()
        << " from region" << regionI
        << " cells back to base mesh." << endl;
    cellProcAddressing.write();

    labelIOList boundaryProcAddressing
    (
        IOobject
        (
            "boundaryRegionAddressing",
            newMesh().facesInstance(),
            polyMesh::meshSubDir,
            newMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        labelList(nNewPatches, -1)
    );
    forAll(oldToNew, i)
    {
        if (!addedPatches.found(i))
        {
            label newI = oldToNew[i];
            if (newI >= 0 && newI < nNewPatches)
            {
                boundaryProcAddressing[oldToNew[i]] = i;
            }
        }
    }
    Info<< "Writing map " << boundaryProcAddressing.name()
        << " from region" << regionI
        << " boundary back to base mesh." << endl;
    boundaryProcAddressing.write();
}


// Create for every region-region interface with non-zero size two patches.
// First one is for minimumregion to maximumregion.
// Note that patches get created in same order on all processors (if parallel)
// since looping over synchronised 'interfaces'.
labelList addRegionPatches
(
    fvMesh& mesh,
    const wordList& regionNames,
    const edgeList& interfaces,
    const List<Pair<word>>& interfaceNames
)
{
    Info<< nl << "Adding patches" << nl << endl;

    labelList interfacePatches(interfaces.size());

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];
        const Pair<word>& names = interfaceNames[interI];

        //Info<< "For interface " << interI
        //    << " between regions " << e
        //    << " trying to add patches " << names << endl;


        mappedWallPolyPatch patch1
        (
            names[0],
            0,                  // overridden
            0,                  // overridden
            0,                  // overridden
            regionNames[e[1]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[1],           // samplePatch
            point::zero,        // offset
            mesh.boundaryMesh()
        );

        interfacePatches[interI] = fvMeshTools::addPatch
        (
            mesh,
            patch1,
            dictionary(),   //optional per field value
            fvPatchFieldBase::calculatedType(),
            true            //validBoundary
        );

        mappedWallPolyPatch patch2
        (
            names[1],
            0,
            0,
            0,
            regionNames[e[0]],  // sampleRegion
            mappedPatchBase::NEARESTPATCHFACE,
            names[0],
            point::zero,        // offset
            mesh.boundaryMesh()
        );
        fvMeshTools::addPatch
        (
            mesh,
            patch2,
            dictionary(),   //optional per field value
            fvPatchFieldBase::calculatedType(),
            true            //validBoundary
        );

        Info<< "For interface between region " << regionNames[e[0]]
            << " and " << regionNames[e[1]] << " added patches" << endl
            << "    " << interfacePatches[interI]
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]].name()
            << endl
            << "    " << interfacePatches[interI]+1
            << "\t" << mesh.boundaryMesh()[interfacePatches[interI]+1].name()
            << endl;
    }
    return interfacePatches;
}


// Find region that covers most of cell zone
label findCorrespondingRegion
(
    const labelList& existingZoneID,    // per cell the (unique) zoneID
    const labelList& cellRegion,
    const label nCellRegions,
    const label zoneI,
    const label minOverlapSize
)
{
    // Per region the number of cells in zoneI
    labelList cellsInZone(nCellRegions, Zero);

    forAll(cellRegion, celli)
    {
        if (existingZoneID[celli] == zoneI)
        {
            cellsInZone[cellRegion[celli]]++;
        }
    }

    Pstream::listCombineReduce(cellsInZone, plusEqOp<label>());

    // Pick region with largest overlap of zoneI
    label regionI = findMax(cellsInZone);


    if (cellsInZone[regionI] < minOverlapSize)
    {
        // Region covers too little of zone. Not good enough.
        regionI = -1;
    }
    else
    {
        // Check that region contains no cells that aren't in cellZone.
        forAll(cellRegion, celli)
        {
            if (cellRegion[celli] == regionI && existingZoneID[celli] != zoneI)
            {
                // celli in regionI but not in zoneI
                regionI = -1;
                break;
            }
        }
        // If one in error, all should be in error. Note that branch gets taken
        // on all procs.
        reduce(regionI, minOp<label>());
    }

    return regionI;
}


void getClusterID
(
    const polyMesh& mesh,
    const cellZoneMesh& cellZones,
    const wordList& clusterNames,
    const labelListList& clusterToZones,
    labelList& clusterID,
    labelList& neiClusterID
)
{
    // Existing zoneID
    clusterID.setSize(mesh.nCells());
    clusterID = -1;

    forAll(clusterToZones, clusterI)
    {
        for (const label zoneI : clusterToZones[clusterI])
        {
            const cellZone& cz = cellZones[zoneI];

            forAll(cz, i)
            {
                label celli = cz[i];
                if (clusterID[celli] == -1)
                {
                    clusterID[celli] = clusterI;
                }
                else
                {
                    FatalErrorInFunction
                        << "Cell " << celli << " with cell centre "
                        << mesh.cellCentres()[celli]
                        << " is multiple zones. This is not allowed." << endl
                        << "It is in zone " << clusterNames[clusterID[celli]]
                        << " and in zone " << clusterNames[clusterI]
                        << exit(FatalError);
                }
            }
        }
    }

    // Neighbour zoneID.
    syncTools::swapBoundaryCellList(mesh, clusterID, neiClusterID);
}


word makeRegionName
(
    const cellZoneMesh& czs,
    const label regioni,
    const labelList& zoneIDs
)
{
    // Synthesise region name. Equals the zone name if cluster consist of only
    // one zone

    if (zoneIDs.empty())
    {
        return word("domain") + Foam::name(regioni);
    }
    else
    {
        // Zone indices are in cellZone order ...
        word regionName(czs[zoneIDs[0]].name());

        // Synthesize name from appended cellZone names
        for (label i = 1; i < zoneIDs.size(); i++)
        {
            regionName += "_" + czs[zoneIDs[i]].name();
        }
        return regionName;
    }
}


void makeClusters
(
    const List<wordRes>& zoneClusters,
    const wordList& zoneClusterNames,
    const cellZoneMesh& cellZones,
    wordList& clusterNames,
    labelListList& clusterToZones,
    labelList& zoneToCluster
)
{
    // Check if there are clustering for zones. If none every zone goes into
    // its own cluster.

    clusterNames.clear();
    clusterToZones.clear();
    zoneToCluster.setSize(cellZones.size());
    zoneToCluster = -1;

    if (zoneClusters.size())
    {
        forAll(zoneClusters, clusteri)
        {
            const labelList zoneIDs(cellZones.indices(zoneClusters[clusteri]));
            UIndirectList<label>(zoneToCluster, zoneIDs) = clusteri;
            clusterNames.append
            (
                zoneClusterNames[clusteri].size()
              ? zoneClusterNames[clusteri]
              : makeRegionName
                (
                    cellZones,
                    clusteri,
                    zoneIDs
                )
            );
            clusterToZones.append(std::move(zoneIDs));
        }

        // Unclustered zone
        forAll(zoneToCluster, zonei)
        {
            if (zoneToCluster[zonei] == -1)
            {
                clusterNames.append(cellZones[zonei].name());
                clusterToZones.append(labelList(1, zonei));
                zoneToCluster[zonei] = clusterToZones.size();
            }
        }
    }
    else
    {
        for (const auto& cellZone : cellZones)
        {
            const label nClusters = clusterToZones.size();
            clusterNames.append(cellZone.name());
            clusterToZones.append(labelList(1, cellZone.index()));
            zoneToCluster[cellZone.index()] = nClusters;
        }
    }
}


void matchRegions
(
    const bool sloppyCellZones,
    const polyMesh& mesh,

    const wordList& clusterNames,
    const labelListList& clusterToZones,
    const labelList& clusterID,

    const label nCellRegions,
    const labelList& cellRegion,

    labelListList& regionToZones,
    wordList& regionNames,
    labelList& zoneToRegion
)
{
    const cellZoneMesh& cellZones = mesh.cellZones();

    regionToZones.setSize(nCellRegions);
    regionToZones = labelList();
    regionNames.setSize(nCellRegions);
    regionNames = word::null;
    zoneToRegion.setSize(cellZones.size(), -1);


    // Sizes per cluster
    labelList clusterSizes(clusterToZones.size(), Zero);
    forAll(clusterToZones, clusterI)
    {
        for (const label zoneI : clusterToZones[clusterI])
        {
            clusterSizes[clusterI] += cellZones[zoneI].size();
        }
        reduce(clusterSizes[clusterI], sumOp<label>());
    }

    if (sloppyCellZones)
    {
        Info<< "Trying to match regions to existing cell zones;"
            << " region can be subset of cell zone." << nl << endl;

        forAll(clusterToZones, clusterI)
        {
            label regionI = findCorrespondingRegion
            (
                clusterID,
                cellRegion,
                nCellRegions,
                clusterI,
                label(0.5*clusterSizes[clusterI]) // minimum overlap
            );

            if (regionI != -1)
            {
                Info<< "Sloppily matched region " << regionI
                    //<< " size " << regionSizes[regionI]
                    << " to cluster " << clusterI
                    << " size " << clusterSizes[clusterI]
                    << endl;
                UIndirectList<label>
                (
                    zoneToRegion,
                    clusterToZones[clusterI]
                ) = regionI;
                regionToZones[regionI] = clusterToZones[clusterI];
                regionNames[regionI] = clusterNames[clusterI];
            }
        }
    }
    else
    {
        Info<< "Trying to match regions to existing cell zones." << nl << endl;

        forAll(clusterToZones, clusterI)
        {
            label regionI = findCorrespondingRegion
            (
                clusterID,
                cellRegion,
                nCellRegions,
                clusterI,
                clusterSizes[clusterI]  // require exact match
            );

            if (regionI != -1)
            {
                UIndirectList<label>
                (
                    zoneToRegion,
                    clusterToZones[clusterI]
                ) = regionI;
                regionToZones[regionI] = clusterToZones[clusterI];
                regionNames[regionI] = clusterNames[clusterI];
            }
        }
    }
    // Allocate region names for unmatched regions.
    forAll(regionNames, regionI)
    {
        if (regionNames[regionI].empty())
        {
            regionNames[regionI] = makeRegionName
            (
                cellZones,
                regionI,
                regionToZones[regionI]
            );
        }
    }
}


void writeCellToRegion(const fvMesh& mesh, const labelList& cellRegion)
{
    // Write to manual decomposition option
    {
        labelIOList cellToRegion
        (
            IOobject
            (
                "cellToRegion",
                mesh.facesInstance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            cellRegion
        );
        cellToRegion.write();

        Info<< "Writing region per cell file (for manual decomposition) to "
            << cellToRegion.objectPath() << nl << endl;
    }
    // Write for postprocessing
    {
        volScalarField cellToRegion
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
            mesh,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        );
        forAll(cellRegion, celli)
        {
            cellToRegion[celli] = cellRegion[celli];
        }
        cellToRegion.write();

        Info<< "Writing region per cell as volScalarField to "
            << cellToRegion.objectPath() << nl << endl;
    }
}




int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Split mesh into multiple regions (detected by walking across faces)"
    );
    #include "addRegionOption.H"
    #include "addOverwriteOption.H"
    argList::addBoolOption
    (
        "cellZones",
        "Additionally split cellZones off into separate regions"
    );
    argList::addBoolOption
    (
        "cellZonesOnly",
        "Use cellZones only to split mesh into regions; do not use walking"
    );
    argList::addOption
    (
        "cellZonesFileOnly",
        "file",
        "Like -cellZonesOnly, but use specified file"
    );
    argList::addOption
    (
        "combineZones",
        "lists of zones",
        "Combine zones in follow-on analysis"
    );
    argList::addOption
    (
        "addZones",
        "lists of zones",
        "Combine zones in follow-on analysis"
    );
    argList::addOption
    (
        "blockedFaces",
        "faceSet",
        "Specify additional region boundaries that walking does not cross"
    );
    argList::addBoolOption
    (
        "makeCellZones",
        "Place cells into cellZones instead of splitting mesh"
    );
    argList::addBoolOption
    (
        "largestOnly",
        "Only write largest region"
    );
    argList::addOption
    (
        "insidePoint",
        "point",
        "Only write region containing point"
    );
    argList::addBoolOption
    (
        "detectOnly",
        "Do not write mesh"
    );
    argList::addBoolOption
    (
        "sloppyCellZones",
        "Try to match heuristically regions to existing cell zones"
    );
    argList::addBoolOption
    (
        "useFaceZones",
        "Use faceZones to patch inter-region faces instead of single patch"
    );
    argList::addBoolOption
    (
        "prefixRegion",
        "Prefix region name to all patches, not just coupling patches"
    );

    argList::noFunctionObjects();  // Never use function objects

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word oldInstance = mesh.pointsInstance();

    word blockedFacesName;
    if (args.readIfPresent("blockedFaces", blockedFacesName))
    {
        Info<< "Reading blocked internal faces from faceSet "
            << blockedFacesName << nl << endl;
    }

    const bool makeCellZones    = args.found("makeCellZones");
    const bool largestOnly      = args.found("largestOnly");
    const bool insidePoint      = args.found("insidePoint");
    const bool useCellZones     = args.found("cellZones");
    const bool useCellZonesOnly = args.found("cellZonesOnly");
    const bool useCellZonesFile = args.found("cellZonesFileOnly");
    const bool combineZones     = args.found("combineZones");
    const bool addZones         = args.found("addZones");
    const bool overwrite        = args.found("overwrite");
    const bool detectOnly       = args.found("detectOnly");
    const bool sloppyCellZones  = args.found("sloppyCellZones");
    const bool useFaceZones     = args.found("useFaceZones");
    const bool prefixRegion     = args.found("prefixRegion");


    if
    (
        (useCellZonesOnly || useCellZonesFile)
     && (useCellZones || blockedFacesName.size())
    )
    {
        FatalErrorInFunction
            << "You cannot specify both -cellZonesOnly or -cellZonesFileOnly"
            << " (which specify complete split)"
            << " in combination with -blockedFaces or -cellZones"
            << " (which imply a split based on topology)"
            << exit(FatalError);
    }


    if (useFaceZones)
    {
        Info<< "Using current faceZones to divide inter-region interfaces"
            << " into multiple patches."
            << nl << endl;
    }
    else
    {
        Info<< "Creating single patch per inter-region interface."
            << nl << endl;
    }



    if (insidePoint && largestOnly)
    {
        FatalErrorInFunction
            << "You cannot specify both -largestOnly"
            << " (keep region with most cells)"
            << " and -insidePoint (keep region containing point)"
            << exit(FatalError);
    }


    // Make sure cellZone names consistent across processors
    mesh.cellZones().checkParallelSync(true);

    List<wordRes> zoneClusters;
    wordList zoneClusterNames;
    if (combineZones)
    {
        if (addZones)
        {
            FatalErrorInFunction
                << "Cannot specify both combineZones and addZones"
                << exit(FatalError);
        }
        zoneClusters = args.get<List<wordRes>>("combineZones");
        zoneClusterNames.setSize(zoneClusters.size());
    }
    else if (addZones)
    {
        zoneClusters = args.get<List<wordRes>>("addZones");
        zoneClusterNames.setSize(zoneClusters.size());
        forAll(zoneClusters, clusteri)
        {
            // Pop of front - is name

            wordRes& wrs = zoneClusters[clusteri];

            zoneClusterNames[clusteri] = wrs[0];

            for (label i = 1; i < wrs.size(); i++)
            {
                wrs[i-1] = wrs[i];
            }
            wrs.setSize(wrs.size()-1);
        }
    }


    // Determine per cell the region it belongs to
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // cellRegion is the labelList with the region per cell.
    labelList cellRegion;
    // Region to zone(s)
    labelListList regionToZones;
    // Name of region
    wordList regionNames;
    // Zone to region
    labelList zoneToRegion;

    label nCellRegions = 0;
    if (useCellZonesOnly)
    {
        Info<< "Using current cellZones to split mesh into regions."
            << " This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        // Collect sets of zones into clusters. If no cluster is just an identity
        // list (cluster 0 is cellZone 0 etc.)
        wordList clusterNames;
        labelListList clusterToZones;
        labelList zoneToCluster;
        makeClusters
        (
            zoneClusters,
            zoneClusterNames,
            mesh.cellZones(),
            clusterNames,
            clusterToZones,
            zoneToCluster
        );

        // Existing clusterID
        labelList clusterID(mesh.nCells(), -1);
        // Neighbour clusterID.
        labelList neiClusterID(mesh.nBoundaryFaces());
        getClusterID
        (
            mesh,
            mesh.cellZones(),
            clusterNames,
            clusterToZones,
            clusterID,
            neiClusterID
        );

        label unzonedCelli = clusterID.find(-1);
        if (unzonedCelli != -1)
        {
            FatalErrorInFunction
                << "For the cellZonesOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCelli
                << " at" << mesh.cellCentres()[unzonedCelli]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = clusterID;
        nCellRegions = gMax(cellRegion)+1;
        zoneToRegion = zoneToCluster;
        regionToZones = clusterToZones;
        regionNames = clusterNames;
    }
    else if (useCellZonesFile)
    {
        const word zoneFile(args["cellZonesFileOnly"]);
        Info<< "Reading split from cellZones file " << zoneFile << endl
            << "This requires all"
            << " cells to be in one and only one cellZone." << nl << endl;

        cellZoneMesh newCellZones
        (
            IOobject
            (
                zoneFile,
                mesh.facesInstance(),
                polyMesh::meshSubDir,
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh
        );

        wordList clusterNames;
        labelListList clusterToZones;
        labelList zoneToCluster;
        makeClusters
        (
            zoneClusters,
            zoneClusterNames,
            newCellZones,
            clusterNames,
            clusterToZones,
            zoneToCluster
        );


        // Existing clusterID
        labelList clusterID(mesh.nCells(), -1);
        // Neighbour clusterID.
        labelList neiClusterID(mesh.nBoundaryFaces());
        getClusterID
        (
            mesh,
            newCellZones,
            clusterNames,
            clusterToZones,
            clusterID,
            neiClusterID
        );


        label unzonedCelli = clusterID.find(-1);
        if (unzonedCelli != -1)
        {
            FatalErrorInFunction
                << "For the cellZonesFileOnly option all cells "
                << "have to be in a cellZone." << endl
                << "Cell " << unzonedCelli
                << " at" << mesh.cellCentres()[unzonedCelli]
                << " is not in a cellZone. There might be more unzoned cells."
                << exit(FatalError);
        }
        cellRegion = clusterID;
        nCellRegions = gMax(cellRegion)+1;
        zoneToRegion = zoneToCluster;
        regionToZones = clusterToZones;
        regionNames = clusterNames;
    }
    else
    {
        // Determine connected regions
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Mark additional faces that are blocked
        boolList blockedFace;

        // Read from faceSet
        if (blockedFacesName.size())
        {
            faceSet blockedFaceSet(mesh, blockedFacesName);
            Info<< "Read "
                << returnReduce(blockedFaceSet.size(), sumOp<label>())
                << " blocked faces from set " << blockedFacesName << nl << endl;

            blockedFace.setSize(mesh.nFaces(), false);

            for (const label facei : blockedFaceSet)
            {
                blockedFace[facei] = true;
            }
        }

        // Collect sets of zones into clusters. If no cluster is just an
        // identity list (cluster 0 is cellZone 0 etc.)
        wordList clusterNames;
        labelListList clusterToZones;
        labelList zoneToCluster;
        makeClusters
        (
            zoneClusters,
            zoneClusterNames,
            mesh.cellZones(),
            clusterNames,
            clusterToZones,
            zoneToCluster
        );

        // Existing clusterID
        labelList clusterID(mesh.nCells(), -1);
        // Neighbour clusterID.
        labelList neiClusterID(mesh.nBoundaryFaces());
        getClusterID
        (
            mesh,
            mesh.cellZones(),
            clusterNames,
            clusterToZones,
            clusterID,
            neiClusterID
        );


        // Imply from differing cellZones
        if (useCellZones)
        {
            blockedFace.setSize(mesh.nFaces(), false);

            for (label facei = 0; facei < mesh.nInternalFaces(); facei++)
            {
                label ownCluster = clusterID[mesh.faceOwner()[facei]];
                label neiCluster = clusterID[mesh.faceNeighbour()[facei]];

                if (ownCluster != neiCluster)
                {
                    blockedFace[facei] = true;
                }
            }

            // Different cellZones on either side of processor patch.
            forAll(neiClusterID, i)
            {
                label facei = i+mesh.nInternalFaces();
                label ownCluster = clusterID[mesh.faceOwner()[facei]];
                label neiCluster = neiClusterID[i];

                if (ownCluster != neiCluster)
                {
                    blockedFace[facei] = true;
                }
            }
        }

        // Do a topological walk to determine regions
        regionSplit regions(mesh, blockedFace);
        nCellRegions = regions.nRegions();
        cellRegion.transfer(regions);

        // Match regions to zones
        matchRegions
        (
            sloppyCellZones,
            mesh,

            // cluster-to-zone and cluster-to-cell addressing
            clusterNames,
            clusterToZones,
            clusterID,

            // cell-to-region addressing
            nCellRegions,
            cellRegion,

            // matched zones
            regionToZones,
            regionNames,
            zoneToRegion
        );

        // Override any default region names if single region selected
        if (largestOnly || insidePoint)
        {
            forAll(regionToZones, regionI)
            {
                if (regionToZones[regionI].empty())
                {
                    if (overwrite)
                    {
                        regionNames[regionI] = polyMesh::defaultRegion;
                    }
                    else if (insidePoint)
                    {
                        regionNames[regionI] = "insidePoint";
                    }
                    else if (largestOnly)
                    {
                        regionNames[regionI] = "largestOnly";
                    }
                }
            }
        }
    }

    Info<< endl << "Number of regions:" << nCellRegions << nl << endl;


    // Write decomposition to file
    writeCellToRegion(mesh, cellRegion);



    // Sizes per region
    // ~~~~~~~~~~~~~~~~

    labelList regionSizes(nCellRegions, Zero);

    forAll(cellRegion, celli)
    {
        regionSizes[cellRegion[celli]]++;
    }
    forAll(regionSizes, regionI)
    {
        reduce(regionSizes[regionI], sumOp<label>());
    }

    Info<< "Region\tCells" << nl
        << "------\t-----" << endl;

    forAll(regionSizes, regionI)
    {
        Info<< regionI << '\t' << regionSizes[regionI] << nl;
    }
    Info<< endl;



    // Print region to zone
    Info<< "Region\tZone\tName" << nl
        << "------\t----\t----" << endl;
    forAll(regionToZones, regionI)
    {
        Info<< regionI << '\t' << flatOutput(regionToZones[regionI])
            << '\t'
            << regionNames[regionI] << nl;
    }
    Info<< endl;



    // Since we're going to mess with patches and zones make sure all
    // is synchronised
    mesh.boundaryMesh().checkParallelSync(true);
    mesh.faceZones().checkParallelSync(true);


    // Interfaces
    // ----------
    // per interface:
    // - the two regions on either side
    // - the name
    // - the (global) size
    edgeList interfaces;
    List<Pair<word>> interfaceNames;
    labelList interfaceSizes;
    // per face the interface
    labelList faceToInterface;

    getInterfaceSizes
    (
        mesh,
        useFaceZones,
        cellRegion,
        regionNames,

        interfaces,
        interfaceNames,
        interfaceSizes,
        faceToInterface
    );

    Info<< "Sizes of interfaces between regions:" << nl << nl
        << "Interface\tRegion\tRegion\tFaces" << nl
        << "---------\t------\t------\t-----" << endl;

    forAll(interfaces, interI)
    {
        const edge& e = interfaces[interI];

        Info<< interI
            << "\t\t" << e[0] << "\t" << e[1]
            << "\t" << interfaceSizes[interI] << nl;
    }
    Info<< endl;


    if (detectOnly)
    {
        Info<< "End\n" << endl;

        return 0;
    }


    // Read objects in time directory
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    IOobjectList objects(mesh, runTime.timeName());

    // Read vol fields.

    PtrList<volScalarField> vsFlds;
    ReadFields(mesh, objects, vsFlds);

    PtrList<volVectorField> vvFlds;
    ReadFields(mesh, objects, vvFlds);

    PtrList<volSphericalTensorField> vstFlds;
    ReadFields(mesh, objects, vstFlds);

    PtrList<volSymmTensorField> vsymtFlds;
    ReadFields(mesh, objects, vsymtFlds);

    PtrList<volTensorField> vtFlds;
    ReadFields(mesh, objects, vtFlds);

    // Read surface fields.

    PtrList<surfaceScalarField> ssFlds;
    ReadFields(mesh, objects, ssFlds);

    PtrList<surfaceVectorField> svFlds;
    ReadFields(mesh, objects, svFlds);

    PtrList<surfaceSphericalTensorField> sstFlds;
    ReadFields(mesh, objects, sstFlds);

    PtrList<surfaceSymmTensorField> ssymtFlds;
    ReadFields(mesh, objects, ssymtFlds);

    PtrList<surfaceTensorField> stFlds;
    ReadFields(mesh, objects, stFlds);

    Info<< endl;


    // Remove any demand-driven fields ('S', 'V' etc)
    mesh.clearOut();


    if (nCellRegions == 1)
    {
        Info<< "Only one region. Doing nothing." << endl;
    }
    else if (makeCellZones)
    {
        Info<< "Putting cells into cellZones instead of splitting mesh."
            << endl;

        // Check if region overlaps with existing zone. If so keep.

        for (label regionI = 0; regionI < nCellRegions; regionI++)
        {
            const labelList& zones = regionToZones[regionI];

            if (zones.size() == 1 && zones[0] != -1)
            {
                // Existing zone
                const label zoneI = zones[0];
                Info<< "    Region " << regionI << " : corresponds to existing"
                    << " cellZone "
                    << zoneI << ' ' << mesh.cellZones()[zoneI].name() << endl;
            }
            else
            {
                // Create new cellZone.
                const labelList regionCells(findIndices(cellRegion, regionI));

                const word zoneName = "region" + Foam::name(regionI);

                label zoneI = mesh.cellZones().findZoneID(zoneName);

                if (zoneI == -1)
                {
                    zoneI = mesh.cellZones().size();
                    mesh.cellZones().setSize(zoneI+1);
                    mesh.cellZones().set
                    (
                        zoneI,
                        new cellZone
                        (
                            zoneName,           //name
                            std::move(regionCells), //addressing
                            zoneI,              //index
                            mesh.cellZones()    //cellZoneMesh
                        )
                    );
                }
                else
                {
                    mesh.cellZones()[zoneI].clearAddressing();
                    mesh.cellZones()[zoneI] = regionCells;
                }
                Info<< "    Region " << regionI << " : created new cellZone "
                    << zoneI << ' ' << mesh.cellZones()[zoneI].name() << endl;
            }
        }
        mesh.cellZones().writeOpt(IOobject::AUTO_WRITE);

        if (!overwrite)
        {
            ++runTime;
            mesh.setInstance(runTime.timeName());
        }
        else
        {
            mesh.setInstance(oldInstance);
        }

        Info<< "Writing cellZones as new mesh to time " << runTime.timeName()
            << nl << endl;

        mesh.write();


        // Write cellSets for convenience
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        Info<< "Writing cellSets corresponding to cellZones." << nl << endl;

        for (const auto& cz : mesh.cellZones())
        {
            cellSet(mesh, cz.name(), cz).write();
        }
    }
    else
    {
        // Add patches for interfaces
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~

        // Add all possible patches. Empty ones get filtered later on.
        Info<< nl << "Adding patches" << nl << endl;

        labelList interfacePatches
        (
            addRegionPatches
            (
                mesh,
                regionNames,
                interfaces,
                interfaceNames
            )
        );


        if (!overwrite)
        {
            ++runTime;
        }


        // Create regions
        // ~~~~~~~~~~~~~~

        if (insidePoint)
        {
            const point insidePoint = args.get<point>("insidePoint");

            label regionI = -1;

            (void)mesh.tetBasePtIs();

            label celli = mesh.findCell(insidePoint);

            Info<< nl << "Found point " << insidePoint << " in cell " << celli
                << endl;

            if (celli != -1)
            {
                regionI = cellRegion[celli];
            }

            reduce(regionI, maxOp<label>());

            Info<< nl
                << "Subsetting region " << regionI
                << " containing point " << insidePoint << endl;

            if (regionI == -1)
            {
                FatalErrorInFunction
                    << "Point " << insidePoint
                    << " is not inside the mesh." << nl
                    << "Bounding box of the mesh:" << mesh.bounds()
                    << exit(FatalError);
            }

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                prefixRegion,
                faceToInterface,
                interfacePatches,
                regionI,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else if (largestOnly)
        {
            label regionI = findMax(regionSizes);

            Info<< nl
                << "Subsetting region " << regionI
                << " of size " << regionSizes[regionI]
                << " as named region " << regionNames[regionI] << endl;

            createAndWriteRegion
            (
                mesh,
                cellRegion,
                regionNames,
                prefixRegion,
                faceToInterface,
                interfacePatches,
                regionI,
                (overwrite ? oldInstance : runTime.timeName())
            );
        }
        else
        {
            // Split all
            for (label regionI = 0; regionI < nCellRegions; regionI++)
            {
                Info<< nl
                    << "Region " << regionI << nl
                    << "-------- " << endl;

                createAndWriteRegion
                (
                    mesh,
                    cellRegion,
                    regionNames,
                    prefixRegion,
                    faceToInterface,
                    interfacePatches,
                    regionI,
                    (overwrite ? oldInstance : runTime.timeName())
                );
            }
        }
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
