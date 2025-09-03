/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "fluxSummary.H"
#include "surfaceFields.H"
#include "polySurfaceFields.H"
#include "dictionary.H"
#include "Time.H"
#include "syncTools.H"
#include "meshTools.H"
#include "PatchEdgeFaceWave.H"
#include "edgeTopoDistanceData.H"
#include "globalIndex.H"
#include "OBJstream.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fluxSummary, 0);
    addToRunTimeSelectionTable(functionObject, fluxSummary, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::fluxSummary::modeType
>
Foam::functionObjects::fluxSummary::modeTypeNames_
({
    { modeType::mdFaceZone , "faceZone" },
    { modeType::mdFaceZoneAndDirection, "faceZoneAndDirection" },
    { modeType::mdCellZoneAndDirection, "cellZoneAndDirection" },
    { modeType::mdSurface, "functionObjectSurface" },
    { modeType::mdSurface, "surface" },
    { modeType::mdSurfaceAndDirection, "surfaceAndDirection" },
});


// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * * //

bool Foam::functionObjects::fluxSummary::isSurfaceMode() const
{
    return (mdSurface == mode_ || mdSurfaceAndDirection == mode_);
}


Foam::word Foam::functionObjects::fluxSummary::checkFlowType
(
    const dimensionSet& fieldDims,
    const word& fieldName
) const
{
    // Surfaces are multiplied by their area, so account for that
    // in the dimension checking
    const dimensionSet dims =
        (fieldDims * (isSurfaceMode() ? dimTime*dimArea : dimTime));

    if (dims == dimVolume)
    {
        return "volumetric";
    }
    else if (dims == dimMass)
    {
        return "mass";
    }

    FatalErrorInFunction
        << "Unsupported flux field " << fieldName << " with dimensions "
        << fieldDims
        << ".  Expected either mass flow or volumetric flow rate."
        << abort(FatalError);

    return word::null;
}


void Foam::functionObjects::fluxSummary::initialiseSurface
(
    const word& surfName,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<boolList>& faceFlip
) const
{
    const polySurface* surfptr =
        storedObjects().cfindObject<polySurface>(surfName);

    if (!surfptr)
    {
        FatalErrorInFunction
            << "Unable to find surface " << surfName
            << ".  Valid surfaces: "
            << storedObjects().sortedNames<polySurface>() << nl
            << exit(FatalError);
    }

    names.append(surfName);
    directions.append(Zero);        // Dummy value
    faceFlip.append(boolList());    // No flip-map
}


void Foam::functionObjects::fluxSummary::initialiseSurfaceAndDirection
(
    const word& surfName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<boolList>& faceFlip
) const
{
    const polySurface* surfptr =
        storedObjects().cfindObject<polySurface>(surfName);

    if (!surfptr)
    {
        FatalErrorInFunction
            << "Unable to find surface " << surfName
            << ".  Valid surfaces: "
            << storedObjects().sortedNames<polySurface>() << nl
            << exit(FatalError);
    }

    const auto& s = *surfptr;
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    names.append(surfName);
    directions.append(refDir);
    faceFlip.append(boolList());    // No flip-map

    boolList& flips = faceFlip[faceFlip.size()-1];
    flips.setSize(s.size(), false);

    forAll(s, i)
    {
        // Orientation set by comparison with reference direction
        const vector& n = s.faceNormals()[i];

        if ((n & refDir) > tolerance_)
        {
            flips[i] = false;
        }
        else
        {
            flips[i] = true;
        }
    }
}


void Foam::functionObjects::fluxSummary::initialiseFaceZone
(
    const word& faceZoneName,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<labelList>& faceID,
    DynamicList<labelList>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    label zonei = mesh_.faceZones().findZoneID(faceZoneName);
    if (zonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid zones: "
            << mesh_.faceZones().sortedNames() << nl
            << exit(FatalError);
    }
    const faceZone& fZone = mesh_.faceZones()[zonei];

    names.append(faceZoneName);
    directions.append(Zero); // dummy value

    labelList faceIds(fZone.size());
    labelList facePatchIds(fZone.size());
    boolList  faceFlips(fZone.size());

    // Total number of faces selected
    label numFaces = 0;

    forAll(fZone, i)
    {
        const label meshFacei = fZone[i];
        const bool isFlip = fZone.flipMap()[i];

        // Internal faces
        label faceId = meshFacei;
        label facePatchId = -1;

        // Boundary faces
        if (!mesh_.isInternalFace(meshFacei))
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(meshFacei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];

            if (isA<emptyPolyPatch>(pp))
            {
                continue;  // Ignore empty patch
            }

            const auto* cpp = isA<coupledPolyPatch>(pp);

            if (cpp && !cpp->owner())
            {
                continue;  // Ignore neighbour side
            }

            faceId = pp.whichFace(meshFacei);
        }

        if (faceId >= 0)
        {
            // Orientation set by faceZone flip map
            faceIds[numFaces] = faceId;
            facePatchIds[numFaces] = facePatchId;
            faceFlips[numFaces] = isFlip;

            ++numFaces;
        }
    }

    // Shrink to size used
    faceIds.resize(numFaces);
    facePatchIds.resize(numFaces);
    faceFlips.resize(numFaces);

    faceID.append(std::move(faceIds));
    facePatchID.append(std::move(facePatchIds));
    faceFlip.append(std::move(faceFlips));
}


void Foam::functionObjects::fluxSummary::initialiseFaceZoneAndDirection
(
    const word& faceZoneName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<labelList>& faceID,
    DynamicList<labelList>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    label zonei = mesh_.faceZones().findZoneID(faceZoneName);
    if (zonei == -1)
    {
         FatalErrorInFunction
            << "Unable to find faceZone " << faceZoneName
            << ".  Valid zones: "
            << mesh_.faceZones().sortedNames() << nl
            << exit(FatalError);
    }
    const faceZone& fZone = mesh_.faceZones()[zonei];

    names.append(faceZoneName);
    directions.append(refDir);

    labelList faceIds(fZone.size());
    labelList facePatchIds(fZone.size());
    boolList  faceFlips(fZone.size());

    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    vector n(Zero);

    // Total number of faces selected
    label numFaces = 0;

    forAll(fZone, i)
    {
        const label meshFacei = fZone[i];

        // Internal faces
        label faceId = meshFacei;
        label facePatchId = -1;

        // Boundary faces
        if (!mesh_.isInternalFace(meshFacei))
        {
            facePatchId = mesh_.boundaryMesh().whichPatch(meshFacei);
            const polyPatch& pp = mesh_.boundaryMesh()[facePatchId];

            if (isA<emptyPolyPatch>(pp))
            {
                continue;  // Ignore empty patch
            }

            const auto* cpp = isA<coupledPolyPatch>(pp);

            if (cpp && !cpp->owner())
            {
                continue;  // Ignore neighbour side
            }

            faceId = pp.whichFace(meshFacei);
        }

        if (faceId >= 0)
        {
            // Orientation set by comparison with reference direction
            if (facePatchId != -1)
            {
                n = Sf.boundaryField()[facePatchId][faceId]
                   /(magSf.boundaryField()[facePatchId][faceId] + ROOTVSMALL);
            }
            else
            {
                n = Sf[faceId]/(magSf[faceId] + ROOTVSMALL);
            }

            if ((n & refDir) > tolerance_)
            {
                faceFlips[numFaces] = false;
            }
            else
            {
                faceFlips[numFaces] = true;
            }

            faceIds[numFaces] = faceId;
            facePatchIds[numFaces] = facePatchId;

            ++numFaces;
        }
    }

    // Shrink to size used
    faceIds.resize(numFaces);
    facePatchIds.resize(numFaces);
    faceFlips.resize(numFaces);

    faceID.append(std::move(faceIds));
    facePatchID.append(std::move(facePatchIds));
    faceFlip.append(std::move(faceFlips));
}


void Foam::functionObjects::fluxSummary::initialiseCellZoneAndDirection
(
    const word& cellZoneName,
    const vector& dir,
    DynamicList<word>& names,
    DynamicList<vector>& directions,
    DynamicList<labelList>& faceID,
    DynamicList<labelList>& facePatchID,
    DynamicList<boolList>& faceFlip
) const
{
    const vector refDir = dir/(mag(dir) + ROOTVSMALL);

    const label cellZonei = mesh_.cellZones().findZoneID(cellZoneName);
    if (cellZonei == -1)
    {
        FatalErrorInFunction
            << "Unable to find cellZone " << cellZoneName
            << ".  Valid zones: "
            << mesh_.cellZones().sortedNames() << nl
            << exit(FatalError);
    }

    const label nInternalFaces = mesh_.nInternalFaces();
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    labelList cellAddr(mesh_.nCells(), -1);
    const labelList& cellIDs = mesh_.cellZones()[cellZonei];
    labelUIndList(cellAddr, cellIDs) = identity(cellIDs.size());
    labelList nbrFaceCellAddr(mesh_.nBoundaryFaces(), -1);

    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        if (pp.coupled())
        {
            forAll(pp, i)
            {
                label facei = pp.start() + i;
                label nbrFacei = facei - nInternalFaces;
                label own = mesh_.faceOwner()[facei];
                nbrFaceCellAddr[nbrFacei] = cellAddr[own];
            }
        }
    }

    // Correct boundary values for parallel running
    syncTools::swapBoundaryFaceList(mesh_, nbrFaceCellAddr);

    // Collect faces
    DynamicList<label> faceIDs(floor(0.1*mesh_.nFaces()));
    DynamicList<label> facePatchIDs(faceIDs.size());
    DynamicList<label> faceLocalPatchIDs(faceIDs.size());
    DynamicList<bool>  flips(faceIDs.size());

    // Internal faces
    for (label facei = 0; facei < nInternalFaces; facei++)
    {
        const label own = cellAddr[mesh_.faceOwner()[facei]];
        const label nbr = cellAddr[mesh_.faceNeighbour()[facei]];

        if (((own != -1) && (nbr == -1)) || ((own == -1) && (nbr != -1)))
        {
            vector n = mesh_.faces()[facei].unitNormal(mesh_.points());

            if ((n & refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                flips.append(false);
            }
            else if ((n & -refDir) > tolerance_)
            {
                faceIDs.append(facei);
                faceLocalPatchIDs.append(facei);
                facePatchIDs.append(-1);
                flips.append(true);
            }
        }
    }

    // Loop over boundary faces
    forAll(pbm, patchi)
    {
        const polyPatch& pp = pbm[patchi];

        forAll(pp, localFacei)
        {
            const label facei = pp.start() + localFacei;
            const label own = cellAddr[mesh_.faceOwner()[facei]];
            const label nbr = nbrFaceCellAddr[facei - nInternalFaces];

            if ((own != -1) && (nbr == -1))
            {
                vector n = mesh_.faces()[facei].unitNormal(mesh_.points());

                if ((n & refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    flips.append(false);
                }
                else if ((n & -refDir) > tolerance_)
                {
                    faceIDs.append(facei);
                    faceLocalPatchIDs.append(localFacei);
                    facePatchIDs.append(patchi);
                    flips.append(true);
                }
            }
        }
    }

    // Convert into primitivePatch for convenience
    indirectPrimitivePatch patch
    (
        IndirectList<face>(mesh_.faces(), faceIDs),
        mesh_.points()
    );

    if (debug)
    {
        OBJstream os(mesh_.time().path()/"patch.obj");
        faceList faces(patch);
        os.write(faces, mesh_.points(), false);
    }


    // Data on all edges and faces
    List<edgeTopoDistanceData<label>> allEdgeInfo(patch.nEdges());
    List<edgeTopoDistanceData<label>> allFaceInfo(patch.size());

    bool search = true;

    DebugInfo
        << "initialiseCellZoneAndDirection: "
        << "Starting walk to split patch into faceZones"
        << endl;

    const globalIndex globalFaces(patch.size());

    label oldFaceID = 0;
    label regioni = 0;

    // Dummy tracking data
    bool dummyData{false};

    while (search)
    {
        DynamicList<label> changedEdges;
        DynamicList<edgeTopoDistanceData<label>> changedInfo;

        label seedFacei = labelMax;
        for (; oldFaceID < patch.size(); oldFaceID++)
        {
            if (!allFaceInfo[oldFaceID].valid<bool>(dummyData))
            {
                seedFacei = globalFaces.toGlobal(oldFaceID);
                break;
            }
        }
        reduce(seedFacei, minOp<label>());

        if (seedFacei == labelMax)
        {
            break;
        }

        if (globalFaces.isLocal(seedFacei))
        {
            const label localFacei = globalFaces.toLocal(seedFacei);
            const labelList& fEdges = patch.faceEdges()[localFacei];

            for (const label edgei : fEdges)
            {
                if (allEdgeInfo[edgei].valid<bool>(dummyData))
                {
                    WarningInFunction
                        << "Problem in edge face wave: attempted to assign a "
                        << "value to an edge that has already been visited. "
                        << "Edge info: " << allEdgeInfo[edgei]
                        << endl;
                }

                changedEdges.append(edgei);
                changedInfo.append
                (
                    edgeTopoDistanceData<label>
                    (
                        0,          // distance
                        regioni
                    )
                );
            }
        }


        PatchEdgeFaceWave
        <
            indirectPrimitivePatch,
            edgeTopoDistanceData<label>
        > calc
        (
            mesh_,
            patch,
            changedEdges,
            changedInfo,
            allEdgeInfo,
            allFaceInfo,
            returnReduce(patch.nEdges(), sumOp<label>())
        );

        if (debug)
        {
            label nCells = 0;
            forAll(allFaceInfo, facei)
            {
                if (allFaceInfo[facei].data() == regioni)
                {
                    nCells++;
                }
            }

            Info<< "*** region:" << regioni
                << "  found:" << returnReduce(nCells, sumOp<label>())
                << " faces" << endl;
        }

        regioni++;
    }

    // Collect the data per region
    const label nRegion = regioni;

    #ifdef FULLDEBUG
    // May wish to have enabled always
    if (nRegion == 0)
    {
         FatalErrorInFunction
            << "Region split failed" << nl
            << exit(FatalError);
    }
    #endif

    List<DynamicList<label>> regionFaceIDs(nRegion);
    List<DynamicList<label>> regionFacePatchIDs(nRegion);
    List<DynamicList<bool>>  regionFaceFlips(nRegion);

    forAll(allFaceInfo, facei)
    {
        regioni = allFaceInfo[facei].data();

        regionFaceIDs[regioni].append(faceLocalPatchIDs[facei]);
        regionFacePatchIDs[regioni].append(facePatchIDs[facei]);
        regionFaceFlips[regioni].append(flips[facei]);
    }

    // Transfer to persistent storage
    forAll(regionFaceIDs, regioni)
    {
        const word zoneName = cellZoneName + ":faceZone" + Foam::name(regioni);
        names.append(zoneName);
        directions.append(refDir);
        faceID.append(regionFaceIDs[regioni]);
        facePatchID.append(regionFacePatchIDs[regioni]);
        faceFlip.append(regionFaceFlips[regioni]);

        // Write OBj of faces to file
        if (debug)
        {
            OBJstream os(mesh_.time().path()/zoneName + ".obj");
            faceList faces(mesh_.faces(), regionFaceIDs[regioni]);
            os.write(faces, mesh_.points(), false);
        }
    }

    if (log)
    {
        Info<< type() << " " << name() << " output:" << nl
            << "    Created " << faceID.size()
            << " separate face zones from cell zone " << cellZoneName << nl;

        forAll(names, i)
        {
            label nFaces = returnReduce(faceID[i].size(), sumOp<label>());
            Info<< "    " << names[i] << ": "
                << nFaces << " faces" << nl;
        }

        Info<< endl;
    }
}


Foam::scalar Foam::functionObjects::fluxSummary::totalArea
(
    const label idx
) const
{
    scalar sumMagSf = 0;

    if (isSurfaceMode())
    {
        const polySurface& s =
            storedObjects().lookupObject<polySurface>(zoneNames_[idx]);

        sumMagSf = sum(s.magSf());
    }
    else
    {
        const surfaceScalarField& magSf = mesh_.magSf();

        const labelList& faceIDs = faceID_[idx];
        const labelList& facePatchIDs = facePatchID_[idx];

        forAll(faceIDs, i)
        {
            label facei = faceIDs[i];

            if (facePatchIDs[i] == -1)
            {
                sumMagSf += magSf[facei];
            }
            else
            {
                label patchi = facePatchIDs[i];
                sumMagSf += magSf.boundaryField()[patchi][facei];
            }
        }
    }

    return returnReduce(sumMagSf, sumOp<scalar>());
}


bool Foam::functionObjects::fluxSummary::surfaceModeWrite()
{
    for (const word& surfName : zoneNames_)
    {
        const polySurface& s =
            storedObjects().lookupObject<polySurface>(surfName);

        const auto& phi = s.lookupObject<polySurfaceVectorField>(phiName_);

        Log << type() << ' ' << name() << ' '
            << checkFlowType(phi.dimensions(), phi.name()) << " write:" << nl;
    }


    forAll(zoneNames_, surfi)
    {
        const polySurface& s =
            storedObjects().lookupObject<polySurface>(zoneNames_[surfi]);

        const auto& phi = s.lookupObject<polySurfaceVectorField>(phiName_);

        checkFlowType(phi.dimensions(), phi.name());

        const boolList& flips = faceFlip_[surfi];

        scalar phiPos(0);
        scalar phiNeg(0);

        tmp<scalarField> tphis = phi & s.Sf();
        const scalarField& phis = tphis();

        forAll(s, i)
        {
            scalar phif = phis[i];
            if (flips[i])
            {
                phif *= -1;
            }

            if (phif > 0)
            {
                phiPos += phif;
            }
            else
            {
                phiNeg += phif;
            }
        }

        reduce(phiPos, sumOp<scalar>());
        reduce(phiNeg, sumOp<scalar>());

        phiPos *= scaleFactor_;
        phiNeg *= scaleFactor_;

        scalar netFlux = phiPos + phiNeg;
        scalar absoluteFlux = phiPos - phiNeg;

        Log << "    surface " << zoneNames_[surfi] << ':' << nl
            << "        positive : " << phiPos << nl
            << "        negative : " << phiNeg << nl
            << "        net      : " << netFlux << nl
            << "        absolute : " << absoluteFlux
            << nl << endl;

        if (writeToFile())
        {
            filePtrs_[surfi]
                << time_.value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::fluxSummary::update()
{
    if (!needsUpdate_)
    {
        return false;
    }

    // Initialise with capacity == number of input names
    DynamicList<word> faceZoneName(zoneNames_.size());
    DynamicList<vector>    refDir(faceZoneName.capacity());
    DynamicList<labelList> faceID(faceZoneName.capacity());
    DynamicList<labelList> facePatchID(faceZoneName.capacity());
    DynamicList<boolList>  faceFlips(faceZoneName.capacity());

    switch (mode_)
    {
        case mdFaceZone:
        {
            forAll(zoneNames_, zonei)
            {
                initialiseFaceZone
                (
                    zoneNames_[zonei],
                    faceZoneName,
                    refDir, // fill with dummy value
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdFaceZoneAndDirection:
        {
            forAll(zoneNames_, zonei)
            {
                initialiseFaceZoneAndDirection
                (
                    zoneNames_[zonei],
                    zoneDirections_[zonei],
                    faceZoneName,
                    refDir,
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdCellZoneAndDirection:
        {
            forAll(zoneNames_, zonei)
            {
                initialiseCellZoneAndDirection
                (
                    zoneNames_[zonei],
                    zoneDirections_[zonei],
                    faceZoneName,
                    refDir,
                    faceID,
                    facePatchID,
                    faceFlips
                );
            }
            break;
        }
        case mdSurface:
        {
            forAll(zoneNames_, zonei)
            {
                initialiseSurface
                (
                    zoneNames_[zonei],
                    faceZoneName,
                    refDir,
                    faceFlips
                );
            }
            break;
        }
        case mdSurfaceAndDirection:
        {
            forAll(zoneNames_, zonei)
            {
                initialiseSurfaceAndDirection
                (
                    zoneNames_[zonei],
                    zoneDirections_[zonei],
                    faceZoneName,
                    refDir,
                    faceFlips
                );
            }
            break;
        }

        // Compiler warning if we forgot an enumeration
    }

    zoneNames_.transfer(faceZoneName);
    faceID_.transfer(faceID);
    facePatchID_.transfer(facePatchID);
    faceFlip_.transfer(faceFlips);

    Info<< type() << " " << name() << " output:" << nl;

    // Calculate and report areas
    List<scalar> areas(zoneNames_.size());
    forAll(zoneNames_, zonei)
    {
        const word& zoneName = zoneNames_[zonei];
        areas[zonei] = totalArea(zonei);

        if (isSurfaceMode())
        {
            Info<< "    Surface: " << zoneName
                << ", area: " << areas[zonei] << nl;
        }
        else
        {
            Info<< "    Zone: " << zoneName
                << ", area: " << areas[zonei] << nl;
        }
    }
    Info<< endl;

    if (writeToFile())
    {
        filePtrs_.resize(zoneNames_.size());

        forAll(filePtrs_, zonei)
        {
            const word& zoneName = zoneNames_[zonei];
            filePtrs_.set(zonei, newFileAtStartTime(zoneName));
            writeFileHeader
            (
                zoneName,
                areas[zonei],
                refDir[zonei],
                filePtrs_[zonei]
            );
        }
    }

    Info<< endl;

    needsUpdate_ = false;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fluxSummary::fluxSummary
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name),
    needsUpdate_(true),
    mode_(mdFaceZone),
    scaleFactor_(1),
    phiName_("phi"),
    zoneNames_(),
    zoneDirections_(),
    faceID_(),
    facePatchID_(),
    faceFlip_(),
    filePtrs_(),
    tolerance_(0.8)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fluxSummary::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    needsUpdate_ = true;
    mode_ = modeTypeNames_.get("mode", dict);
    phiName_ = dict.getOrDefault<word>("phi", "phi");
    scaleFactor_ = dict.getOrDefault<scalar>("scaleFactor", 1);
    tolerance_ = dict.getOrDefault<scalar>("tolerance", 0.8);

    zoneNames_.clear();
    zoneDirections_.clear();

    List<Tuple2<word, vector>> nameAndDirection;

    switch (mode_)
    {
        case mdFaceZone:
        {
            dict.readEntry("faceZones", zoneNames_);
            break;
        }
        case mdFaceZoneAndDirection:
        {
            dict.readEntry("faceZoneAndDirection", nameAndDirection);
            break;
        }
        case mdCellZoneAndDirection:
        {
            dict.readEntry("cellZoneAndDirection", nameAndDirection);
            break;
        }
        case mdSurface:
        {
            dict.readEntry("surfaces", zoneNames_);
            break;
        }
        case mdSurfaceAndDirection:
        {
            dict.readEntry("surfaceAndDirection", nameAndDirection);
            break;
        }
        default:
        {
            FatalIOErrorInFunction(dict)
                << "unhandled enumeration " << modeTypeNames_[mode_]
                << abort(FatalIOError);
        }
    }


    // Split name/vector into separate lists
    if (nameAndDirection.size())
    {
        zoneNames_.resize(nameAndDirection.size());
        zoneDirections_.resize(nameAndDirection.size());

        label zonei = 0;

        for (const Tuple2<word, vector>& nameDirn : nameAndDirection)
        {
            zoneNames_[zonei] = nameDirn.first();
            zoneDirections_[zonei] = nameDirn.second();
            ++zonei;
        }

        nameAndDirection.clear();
    }


    Info<< type() << ' ' << name() << " ("
        << modeTypeNames_[mode_] << ") with selection:\n    "
        << flatOutput(zoneNames_) << endl;

    return !zoneNames_.empty();
}


void Foam::functionObjects::fluxSummary::writeFileHeader
(
    const word& zoneName,
    const scalar area,
    const vector& refDir,
    Ostream& os
) const
{
    writeHeader(os, "Flux summary");
    if (isSurfaceMode())
    {
        writeHeaderValue(os, "Surface", zoneName);
    }
    else
    {
        writeHeaderValue(os, "Face zone", zoneName);
    }
    writeHeaderValue(os, "Total area", area);

    switch (mode_)
    {
        case mdFaceZoneAndDirection:
        case mdCellZoneAndDirection:
        case mdSurfaceAndDirection:
        {
            writeHeaderValue(os, "Reference direction", refDir);
            break;
        }
        default:
        {}
    }

    writeHeaderValue(os, "Scale factor", scaleFactor_);

    writeCommented(os, "Time");
    os  << tab << "positive"
        << tab << "negative"
        << tab << "net"
        << tab << "absolute"
        << endl;
}


bool Foam::functionObjects::fluxSummary::execute()
{
    return true;
}


bool Foam::functionObjects::fluxSummary::write()
{
    update();

    if (isSurfaceMode())
    {
        return surfaceModeWrite();
    }

    const surfaceScalarField& phi = lookupObject<surfaceScalarField>(phiName_);

    Log << type() << ' ' << name() << ' '
        << checkFlowType(phi.dimensions(), phi.name()) << " write:" << nl;

    forAll(zoneNames_, zonei)
    {
        const labelList& faceID = faceID_[zonei];
        const labelList& facePatchID = facePatchID_[zonei];
        const boolList&  faceFlips = faceFlip_[zonei];

        scalar phiPos(0);
        scalar phiNeg(0);
        scalar phif(0);

        forAll(faceID, i)
        {
            label facei = faceID[i];
            label patchi = facePatchID[i];

            if (patchi != -1)
            {
                phif = phi.boundaryField()[patchi][facei];
            }
            else
            {
                phif = phi[facei];
            }

            if (faceFlips[i])
            {
                phif *= -1;
            }

            if (phif > 0)
            {
                phiPos += phif;
            }
            else
            {
                phiNeg += phif;
            }
        }

        reduce(phiPos, sumOp<scalar>());
        reduce(phiNeg, sumOp<scalar>());

        phiPos *= scaleFactor_;
        phiNeg *= scaleFactor_;

        scalar netFlux = phiPos + phiNeg;
        scalar absoluteFlux = phiPos - phiNeg;

        Log << "    faceZone " << zoneNames_[zonei] << ':' << nl
            << "        positive : " << phiPos << nl
            << "        negative : " << phiNeg << nl
            << "        net      : " << netFlux << nl
            << "        absolute : " << absoluteFlux
            << nl << endl;

        if (writeToFile())
        {
            filePtrs_[zonei]
                << time_.value() << token::TAB
                << phiPos << token::TAB
                << phiNeg << token::TAB
                << netFlux << token::TAB
                << absoluteFlux
                << endl;
        }
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
