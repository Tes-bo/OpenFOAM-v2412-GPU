/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "fvMeshDistribute.H"
#include "fvMeshAdder.H"
#include "fvMeshSubset.H"
#include "faceCoupleInfo.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"
#include "processorCyclicPolyPatch.H"
#include "processorCyclicFvPatchField.H"
#include "polyTopoChange.H"
#include "removeCells.H"
#include "polyModifyFace.H"
#include "polyRemovePoint.H"
#include "mapDistributePolyMesh.H"
#include "surfaceFields.H"
#include "syncTools.H"
#include "CompactListList.H"
#include "fvMeshTools.H"
#include "labelPairHashes.H"
#include "ListOps.H"
#include "globalIndex.H"
#include "cyclicACMIPolyPatch.H"
#include "mappedPatchBase.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fvMeshDistribute, 0);

    //- Less function class that can be used for sorting processor patches
    class lessProcPatches
    {
        const labelList& nbrProc_;
        const labelList& referPatchID_;

    public:

        lessProcPatches(const labelList& nbrProc, const labelList& referPatchID)
        :
            nbrProc_(nbrProc),
            referPatchID_(referPatchID)
        {}

        bool operator()(const label a, const label b)
        {
            if (nbrProc_[a] < nbrProc_[b])
            {
                return true;
            }
            else if (nbrProc_[a] > nbrProc_[b])
            {
                return false;
            }
            else
            {
                // Equal neighbour processor
                return referPatchID_[a] < referPatchID_[b];
            }
        }
    };
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fvMeshDistribute::inplaceRenumberWithFlip
(
    const labelUList& oldToNew,
    const bool oldToNewHasFlip,
    const bool lstHasFlip,
    labelUList& lst
)
{
    if (!lstHasFlip && !oldToNewHasFlip)
    {
        Foam::inplaceRenumber(oldToNew, lst);
    }
    else
    {
        // Either input data or map encodes sign so result encodes sign

        forAll(lst, elemI)
        {
            // Extract old value and sign
            label val = lst[elemI];
            label sign = 1;
            if (lstHasFlip)
            {
                if (val > 0)
                {
                    val = val-1;
                }
                else if (val < 0)
                {
                    val = -val-1;
                    sign = -1;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << val
                        << " at index " << elemI << " out of " << lst.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Lookup new value and possibly change sign
            label newVal = oldToNew[val];

            if (oldToNewHasFlip)
            {
                if (newVal > 0)
                {
                    newVal = newVal-1;
                }
                else if (newVal < 0)
                {
                    newVal = -newVal-1;
                    sign = -sign;
                }
                else
                {
                    FatalErrorInFunction
                        << "Problem : zero value " << newVal
                        << " at index " << elemI << " out of "
                        << oldToNew.size()
                        << " list with flip bit" << exit(FatalError);
                }
            }


            // Encode new value and sign
            lst[elemI] = sign*(newVal+1);
        }
    }
}


Foam::labelList Foam::fvMeshDistribute::select
(
    const bool selectEqual,
    const labelList& values,
    const label value
)
{
    label n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            n++;
        }
    }

    labelList indices(n);
    n = 0;

    forAll(values, i)
    {
        if (selectEqual == (values[i] == value))
        {
            indices[n++] = i;
        }
    }
    return indices;
}


Foam::wordList Foam::fvMeshDistribute::mergeWordList(const wordList& procNames)
{
    List<wordList> allNames(Pstream::nProcs());
    allNames[Pstream::myProcNo()] = procNames;
    Pstream::gatherList(allNames);

    // Assume there are few zone names to merge. Use HashSet otherwise (but
    // maintain ordering)
    DynamicList<word> mergedNames;
    if (Pstream::master())
    {
        mergedNames = procNames;
        for (const wordList& names : allNames)
        {
            for (const word& name : names)
            {
                mergedNames.push_uniq(name);
            }
        }
    }
    Pstream::broadcast(mergedNames);

    return mergedNames;
}


void Foam::fvMeshDistribute::printMeshInfo(const polyMesh& mesh)
{
    Pout<< "Primitives:" << nl
        << "    points       :" << mesh.nPoints() << nl
        << "    bb           :" << boundBox(mesh.points(), false) << nl
        << "    internalFaces:" << mesh.nInternalFaces() << nl
        << "    faces        :" << mesh.nFaces() << nl
        << "    cells        :" << mesh.nCells() << nl;

    Pout<< "Patches:" << endl;
    for (const polyPatch& pp : mesh.boundaryMesh())
    {
        Pout<< "    " << pp.index() << " name:" << pp.name()
            << " size:" << pp.size()
            << " start:" << pp.start()
            << " type:" << pp.type()
            << endl;
    }

    if (mesh.pointZones().empty()) Pout<< "PointZones:" << endl;
    for (const auto& zn : mesh.pointZones())
    {
        Pout<< "    " << zn.index() << " name:" << zn.name()
            << " size:" << zn.size() << endl;
    }

    if (!mesh.faceZones().empty()) Pout<< "FaceZones:" << endl;
    for (const auto& zn : mesh.faceZones())
    {
        Pout<< "    " << zn.index() << " name:" << zn.name()
            << " size:" << zn.size() << endl;
    }

    if (!mesh.cellZones().empty()) Pout<< "CellZones:" << endl;
    for (const auto& zn : mesh.cellZones())
    {
        Pout<< "    " << zn.index() << " name:" << zn.name()
            << " size:" << zn.size() << endl;
    }
}


void Foam::fvMeshDistribute::printCoupleInfo
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc
)
{
    Pout<< nl
        << "Current coupling info:"
        << endl;

    forAll(sourceFace, bFacei)
    {
        label meshFacei = mesh.nInternalFaces() + bFacei;

        Pout<< "    meshFace:" << meshFacei
            << " fc:" << mesh.faceCentres()[meshFacei]
            << " connects to proc:" << sourceProc[bFacei]
            << "/face:" << sourceFace[bFacei]
            << " which will move to proc:" << sourceNewNbrProc[bFacei]
            << endl;
    }
}


Foam::label Foam::fvMeshDistribute::findNonEmptyPatch() const
{
    // Finds (non-empty) patch that exposed internal and proc faces can be
    // put into.
    const polyBoundaryMesh& patches = mesh_.boundaryMesh();


    // Mark 'special' patches : -coupled, -duplicate faces. These probably
    // should not be used to (temporarily) store processor faces ...

    bitSet isCoupledPatch(patches.size());
    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        const auto* cpp = isA<cyclicACMIPolyPatch>(pp);

        if (cpp)
        {
            isCoupledPatch.set(patchi);
            const label dupPatchID = cpp->nonOverlapPatchID();
            if (dupPatchID != -1)
            {
                isCoupledPatch.set(dupPatchID);
            }
        }
        else if (pp.coupled())
        {
            isCoupledPatch.set(patchi);
        }
    }

    label nonEmptyPatchi = -1;

    forAllReverse(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if
        (
           !isA<emptyPolyPatch>(pp)
        && !isCoupledPatch(patchi)
        && !isA<mappedPatchBase>(pp)
        )
        {
            nonEmptyPatchi = patchi;
            break;
        }
    }

    if (nonEmptyPatchi == -1)
    {
        FatalErrorInFunction
            << "Cannot find a patch which is not of type empty, mapped or"
            << " coupled in patches " << patches.names() << endl
            << "There has to be at least one such patch for"
            << " distribution to work" << abort(FatalError);
    }

    if (debug)
    {
        Pout<< "findNonEmptyPatch : using patch " << nonEmptyPatchi
            << " name:" << patches[nonEmptyPatchi].name()
            << " type:" << patches[nonEmptyPatchi].type()
            << " to put exposed faces into." << endl;
    }


    // Do additional test for processor patches intermingled with non-proc
    // patches.
    label procPatchi = -1;

    forAll(patches, patchi)
    {
        if (isA<processorPolyPatch>(patches[patchi]))
        {
            procPatchi = patchi;
        }
        else if (procPatchi != -1)
        {
            FatalErrorInFunction
                << "Processor patches should be at end of patch list."
                << endl
                << "Have processor patch " << procPatchi
                << " followed by non-processor patch " << patchi
                << " in patches " << patches.names()
                << abort(FatalError);
        }
    }

    return nonEmptyPatchi;
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvMeshDistribute::generateTestField
(
    const fvMesh& mesh
)
{
    const vector testNormal = normalised(vector::one);

    tmp<surfaceScalarField> tfld
    (
        new surfaceScalarField
        (
            IOobject
            (
                "myFlux",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        )
    );
    surfaceScalarField& fld = tfld.ref();

    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    forAll(fld, facei)
    {
        fld[facei] = (n[facei] & testNormal);
    }

    surfaceScalarField::Boundary& fluxBf = fld.boundaryFieldRef();
    const surfaceVectorField::Boundary& nBf = n.boundaryField();

    forAll(fluxBf, patchi)
    {
        fvsPatchScalarField& fvp = fluxBf[patchi];

        scalarField newPfld(fvp.size());
        forAll(newPfld, i)
        {
            newPfld[i] = (nBf[patchi][i] & testNormal);
        }
        fvp == newPfld;
    }

    return tfld;
}


void Foam::fvMeshDistribute::testField(const surfaceScalarField& fld)
{
    const fvMesh& mesh = fld.mesh();

    const vector testNormal = normalised(vector::one);

    const surfaceVectorField n(mesh.Sf()/mesh.magSf());

    forAll(fld, facei)
    {
        scalar cos = (n[facei] & testNormal);

        if (mag(cos - fld[facei]) > 1e-6)
        {
            //FatalErrorInFunction
            WarningInFunction
                << "On internal face " << facei << " at "
                << mesh.faceCentres()[facei]
                << " the field value is " << fld[facei]
                << " whereas cos angle of " << testNormal
                << " with mesh normal " << n[facei]
                << " is " << cos
                //<< exit(FatalError);
                << endl;
        }
    }
    forAll(fld.boundaryField(), patchi)
    {
        const fvsPatchScalarField& fvp = fld.boundaryField()[patchi];
        const fvsPatchVectorField& np = n.boundaryField()[patchi];

        forAll(fvp, i)
        {
            scalar cos = (np[i] & testNormal);

            if (mag(cos - fvp[i]) > 1e-6)
            {
                label facei = fvp.patch().start()+i;
                //FatalErrorInFunction
                WarningInFunction
                    << "On face " << facei
                    << " on patch " << fvp.patch().name()
                    << " at " << mesh.faceCentres()[facei]
                    << " the field value is " << fvp[i]
                    << " whereas cos angle of " << testNormal
                    << " with mesh normal " << np[i]
                    << " is " << cos
                    //<< exit(FatalError);
                    << endl;
            }
        }
    }
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::deleteProcPatches
(
    const label destinationPatch
)
{
    // Delete all processor patches. Move any processor faces into the last
    // non-processor patch.

    // New patchID per boundary faces to be repatched. Is -1 (no change)
    // or new patchID
    labelList newPatchID(mesh_.nBoundaryFaces(), -1);

    for (const polyPatch& pp : mesh_.boundaryMesh())
    {
        if (isA<processorPolyPatch>(pp))
        {
            if (debug)
            {
                Pout<< "Moving all faces of patch " << pp.name()
                    << " into patch " << destinationPatch
                    << endl;
            }

            SubList<label>
            (
                newPatchID,
                pp.size(),
                pp.offset()
            ) = destinationPatch;
        }
    }

    // Note: order of boundary faces been kept the same since the
    // destinationPatch is at the end and we have visited the patches in
    // incremental order.
    labelListList dummyFaceMaps;
    autoPtr<mapPolyMesh> map = repatch(newPatchID, dummyFaceMaps);


    // Delete (now empty) processor patches.
    {
        labelList oldToNew(identity(mesh_.boundaryMesh().size()));
        label newi = 0;
        // Non processor patches first
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (!isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newi++;
            }
        }
        label nNonProcPatches = newi;

        // Processor patches as last
        forAll(mesh_.boundaryMesh(), patchi)
        {
            if (isA<processorPolyPatch>(mesh_.boundaryMesh()[patchi]))
            {
                oldToNew[patchi] = newi++;
            }
        }
        fvMeshTools::reorderPatches(mesh_, oldToNew, nNonProcPatches, false);
    }

    return map;
}


Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::repatch
(
    const labelList& newPatchID,         // per boundary face -1 or new patchID
    labelListList& constructFaceMap
)
{
    polyTopoChange meshMod(mesh_);

    forAll(newPatchID, bFacei)
    {
        if (newPatchID[bFacei] != -1)
        {
            label facei = mesh_.nInternalFaces() + bFacei;

            label zoneID = mesh_.faceZones().whichZone(facei);
            bool zoneFlip = false;

            if (zoneID >= 0)
            {
                const faceZone& fZone = mesh_.faceZones()[zoneID];
                zoneFlip = fZone.flipMap()[fZone.whichFace(facei)];
            }

            meshMod.setAction
            (
                polyModifyFace
                (
                    mesh_.faces()[facei],       // modified face
                    facei,                      // label of face
                    mesh_.faceOwner()[facei],   // owner
                    -1,                         // neighbour
                    false,                      // face flip
                    newPatchID[bFacei],         // patch for face
                    false,                      // remove from zone
                    zoneID,                     // zone for face
                    zoneFlip                    // face flip in zone
                )
            );
        }
    }


    // Do mapping of fields from one patchField to the other ourselves since
    // is currently not supported by updateMesh.

    // Store boundary fields (we only do this for surfaceFields)
    PtrList<FieldField<fvsPatchField, scalar>> sFlds;
    saveBoundaryFields<scalar, surfaceMesh>(sFlds);
    PtrList<FieldField<fvsPatchField, vector>> vFlds;
    saveBoundaryFields<vector, surfaceMesh>(vFlds);
    PtrList<FieldField<fvsPatchField, sphericalTensor>> sptFlds;
    saveBoundaryFields<sphericalTensor, surfaceMesh>(sptFlds);
    PtrList<FieldField<fvsPatchField, symmTensor>> sytFlds;
    saveBoundaryFields<symmTensor, surfaceMesh>(sytFlds);
    PtrList<FieldField<fvsPatchField, tensor>> tFlds;
    saveBoundaryFields<tensor, surfaceMesh>(tFlds);

    // Change the mesh (no inflation). Note: parallel comms allowed.
    //
    // NOTE: there is one very particular problem with this ordering.
    // We first create the processor patches and use these to merge out
    // shared points (see mergeSharedPoints below). So temporarily points
    // and edges do not match!

    // TBD: temporarily unset mesh moving to avoid problems in meshflux
    //      mapping. To be investigated.
    const bool oldMoving = mesh_.moving(false);
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false, true);
    mesh_.moving(oldMoving);
    mapPolyMesh& map = *mapPtr;


    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Map patch fields using stored boundary fields. Note: assumes order
    // of fields has not changed in object registry!
    mapBoundaryFields<scalar, surfaceMesh>(map, sFlds);
    mapBoundaryFields<vector, surfaceMesh>(map, vFlds);
    mapBoundaryFields<sphericalTensor, surfaceMesh>(map, sptFlds);
    mapBoundaryFields<symmTensor, surfaceMesh>(map, sytFlds);
    mapBoundaryFields<tensor, surfaceMesh>(map, tFlds);


    // Move mesh (since morphing does not do this)
    if (map.hasMotionPoints())
    {
        mesh_.movePoints(map.preMotionPoints());
    }

    // Adapt constructMaps.

    if (debug)
    {
        label index = map.reverseFaceMap().find(-1);

        if (index != -1)
        {
            FatalErrorInFunction
                << "reverseFaceMap contains -1 at index:"
                << index << endl
                << "This means that the repatch operation was not just"
                << " a shuffle?" << abort(FatalError);
        }
    }

    forAll(constructFaceMap, proci)
    {
        inplaceRenumberWithFlip
        (
            map.reverseFaceMap(),
            false,
            true,
            constructFaceMap[proci]
        );
    }

    return mapPtr;
}


// Detect shared points. Need processor patches to be present.
// Background: when adding bits of mesh one can get points which
// share the same position but are only detectable to be topologically
// the same point when doing parallel analysis. This routine will
// merge those points.
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::mergeSharedPoints
(
    const labelList& pointToGlobalMaster,
    labelListList& constructPointMap
)
{
    // Find out which sets of points get merged and create a map from
    // mesh point to unique point.

    label nShared = 0;
    forAll(pointToGlobalMaster, pointi)
    {
        if (pointToGlobalMaster[pointi] != -1)
        {
            nShared++;
        }
    }

    if (debug)
    {
        Pout<< "mergeSharedPoints : found " << nShared
            << " points on processor boundaries" << nl << endl;
    }

    Map<label> globalMasterToLocalMaster(2*nShared);
    Map<label> pointToMaster(2*nShared);
    label nMatch = 0;

    forAll(pointToGlobalMaster, pointi)
    {
        label globali = pointToGlobalMaster[pointi];
        if (globali != -1)
        {
            const auto iter = globalMasterToLocalMaster.cfind(globali);

            if (iter.good())
            {
                // Matched to existing master
                pointToMaster.insert(pointi, *iter);
                nMatch++;
            }
            else
            {
                // Found first point. Designate as master
                globalMasterToLocalMaster.insert(globali, pointi);
                pointToMaster.insert(pointi, pointi);
            }
        }
    }

    reduce(nMatch, sumOp<label>());

    if (debug)
    {
        Pout<< "mergeSharedPoints : found "
            << nMatch << " mergeable points" << nl << endl;
    }


    if (nMatch == 0)
    {
        return nullptr;
    }


    polyTopoChange meshMod(mesh_);

    fvMeshAdder::mergePoints(mesh_, pointToMaster, meshMod);

    // Change the mesh (no inflation). Note: parallel comms allowed.

    // TBD: temporarily unset mesh moving to avoid problems in meshflux
    //      mapping. To be investigated.
    const bool oldMoving = mesh_.moving(false);
    autoPtr<mapPolyMesh> mapPtr = meshMod.changeMesh(mesh_, false, true);
    mesh_.moving(oldMoving);
    mapPolyMesh& map = *mapPtr;

    // Update fields. No inflation, parallel sync.
    mesh_.updateMesh(map);

    // Adapt constructMaps for merged points.
    forAll(constructPointMap, proci)
    {
        labelList& constructMap = constructPointMap[proci];

        forAll(constructMap, i)
        {
            label oldPointi = constructMap[i];

            label newPointi = map.reversePointMap()[oldPointi];

            if (newPointi < -1)
            {
                constructMap[i] = -newPointi-2;
            }
            else if (newPointi >= 0)
            {
                constructMap[i] = newPointi;
            }
            else
            {
                FatalErrorInFunction
                    << "Problem. oldPointi:" << oldPointi
                    << " newPointi:" << newPointi << abort(FatalError);
            }
        }
    }

    return mapPtr;
}


void Foam::fvMeshDistribute::getCouplingData
(
    const labelList& distribution,
    labelList& sourceFace,
    labelList& sourceProc,
    labelList& sourcePatch,
    labelList& sourceNewNbrProc,
    labelList& sourcePointMaster
) const
{
    // Construct the coupling information for all (boundary) faces and
    // points

    const label nBnd = mesh_.nBoundaryFaces();
    sourceFace.setSize(nBnd);
    sourceProc.setSize(nBnd);
    sourcePatch.setSize(nBnd);
    sourceNewNbrProc.setSize(nBnd);

    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Get neighbouring meshFace labels and new processor of coupled boundaries.
    labelList nbrFaces(nBnd, -1);
    labelList nbrNewNbrProc(nBnd, -1);

    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];

        if (pp.coupled())
        {
            label offset = pp.start() - mesh_.nInternalFaces();

            // Mesh labels of faces on this side
            forAll(pp, i)
            {
                label bndI = offset + i;
                nbrFaces[bndI] = pp.start()+i;
            }

            // Which processor they will end up on
            SubList<label>(nbrNewNbrProc, pp.size(), offset) =
                labelUIndList(distribution, pp.faceCells())();
        }
    }


    // Exchange the boundary data
    syncTools::swapBoundaryFaceList(mesh_, nbrFaces);
    syncTools::swapBoundaryFaceList(mesh_, nbrNewNbrProc);


    forAll(patches, patchi)
    {
        const polyPatch& pp = patches[patchi];
        label offset = pp.start() - mesh_.nInternalFaces();

        if (isA<processorPolyPatch>(pp))
        {
            const processorPolyPatch& procPatch =
                refCast<const processorPolyPatch>(pp);

            // Check which of the two faces we store.

            if (procPatch.owner())
            {
                // Use my local face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                // Use my neighbours face labels
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = procPatch.neighbProcNo();
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }


            label patchi = -1;
            if (isA<processorCyclicPolyPatch>(pp))
            {
                patchi = refCast<const processorCyclicPolyPatch>
                (
                    pp
                ).referPatchID();
            }

            forAll(pp, i)
            {
                label bndI = offset + i;
                sourcePatch[bndI] = patchi;
            }
        }
        else if (isA<cyclicPolyPatch>(pp))
        {
            const cyclicPolyPatch& cpp = refCast<const cyclicPolyPatch>(pp);

            if (cpp.owner())
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = pp.start()+i;
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
            else
            {
                forAll(pp, i)
                {
                    label bndI = offset + i;
                    sourceFace[bndI] = nbrFaces[bndI];
                    sourceProc[bndI] = Pstream::myProcNo();
                    sourcePatch[bndI] = patchi;
                    sourceNewNbrProc[bndI] = nbrNewNbrProc[bndI];
                }
            }
        }
        else
        {
            // Normal (physical) boundary
            forAll(pp, i)
            {
                label bndI = offset + i;
                sourceFace[bndI] = -1;
                sourceProc[bndI] = -1;
                sourcePatch[bndI] = patchi;
                sourceNewNbrProc[bndI] = -1;
            }
        }
    }


    // Collect coupled (collocated) points
    sourcePointMaster.resize_nocopy(mesh_.nPoints());
    sourcePointMaster = -1;
    {
        // Assign global master point
        const globalIndex globalPoints(mesh_.nPoints());

        const globalMeshData& gmd = mesh_.globalData();
        const indirectPrimitivePatch& cpp = gmd.coupledPatch();
        const labelList& meshPoints = cpp.meshPoints();
        const mapDistribute& slavesMap = gmd.globalCoPointSlavesMap();
        const labelListList& slaves = gmd.globalCoPointSlaves();

        labelList elems(slavesMap.constructSize(), -1);
        forAll(meshPoints, pointi)
        {
            const labelList& slots = slaves[pointi];

            if (slots.size())
            {
                // pointi is a master. Assign a unique label.

                label globalPointi = globalPoints.toGlobal(meshPoints[pointi]);
                elems[pointi] = globalPointi;
                forAll(slots, i)
                {
                    label sloti = slots[i];
                    if (sloti >= meshPoints.size())
                    {
                        // Filter out local collocated points. We don't want
                        // to merge these
                        elems[slots[i]] = globalPointi;
                    }
                }
            }
        }

        // Push slave-slot data back to slaves
        slavesMap.reverseDistribute(elems.size(), elems, false);

        // Extract back onto mesh
        forAll(meshPoints, pointi)
        {
            sourcePointMaster[meshPoints[pointi]] = elems[pointi];
        }
    }
}


// Subset the neighbourCell/neighbourProc fields
void Foam::fvMeshDistribute::subsetCouplingData
(
    const fvMesh& mesh,
    const labelList& pointMap,
    const labelList& faceMap,
    const labelList& cellMap,

    const labelList& oldDistribution,
    const labelList& oldFaceOwner,
    const labelList& oldFaceNeighbour,
    const label oldInternalFaces,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,

    labelList& subFace,
    labelList& subProc,
    labelList& subPatch,
    labelList& subNewNbrProc,
    labelList& subPointMaster
)
{
    subFace.setSize(mesh.nBoundaryFaces());
    subProc.setSize(mesh.nBoundaryFaces());
    subPatch.setSize(mesh.nBoundaryFaces());
    subNewNbrProc.setSize(mesh.nBoundaryFaces());

    forAll(subFace, newBFacei)
    {
        label newFacei = newBFacei + mesh.nInternalFaces();

        label oldFacei = faceMap[newFacei];

        // Was oldFacei internal face? If so which side did we get.
        if (oldFacei < oldInternalFaces)
        {
            subFace[newBFacei] = oldFacei;
            subProc[newBFacei] = Pstream::myProcNo();
            subPatch[newBFacei] = -1;

            label oldOwn = oldFaceOwner[oldFacei];
            label oldNei = oldFaceNeighbour[oldFacei];

            if (oldOwn == cellMap[mesh.faceOwner()[newFacei]])
            {
                // We kept the owner side. Where does the neighbour move to?
                subNewNbrProc[newBFacei] = oldDistribution[oldNei];
            }
            else
            {
                // We kept the neighbour side.
                subNewNbrProc[newBFacei] = oldDistribution[oldOwn];
            }
        }
        else
        {
            // Was boundary face. Take over boundary information
            label oldBFacei = oldFacei - oldInternalFaces;

            subFace[newBFacei] = sourceFace[oldBFacei];
            subProc[newBFacei] = sourceProc[oldBFacei];
            subPatch[newBFacei] = sourcePatch[oldBFacei];
            subNewNbrProc[newBFacei] = sourceNewNbrProc[oldBFacei];
        }
    }


    subPointMaster = UIndirectList<label>(sourcePointMaster, pointMap);
}


// Find cells on mesh whose faceID/procID match the neighbour cell/proc of
// domainMesh. Store the matching face.
void Foam::fvMeshDistribute::findCouples
(
    const primitiveMesh& mesh,
    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,

    const label domain,
    const primitiveMesh& domainMesh,
    const labelList& domainFace,
    const labelList& domainProc,
    const labelList& domainPatch,

    labelList& masterCoupledFaces,
    labelList& slaveCoupledFaces
)
{
    // Store domain neighbour as map so we can easily look for pair
    // with same face+proc.
    labelPairLookup map(domainFace.size());

    forAll(domainProc, bFacei)
    {
        if (domainProc[bFacei] != -1 && domainPatch[bFacei] == -1)
        {
            map.insert
            (
                labelPair(domainFace[bFacei], domainProc[bFacei]),
                bFacei
            );
        }
    }


    // Try to match mesh data.

    masterCoupledFaces.setSize(domainFace.size());
    slaveCoupledFaces.setSize(domainFace.size());
    label coupledI = 0;

    forAll(sourceFace, bFacei)
    {
        if (sourceProc[bFacei] != -1 && sourcePatch[bFacei] == -1)
        {
            labelPair myData(sourceFace[bFacei], sourceProc[bFacei]);

            const auto iter = map.cfind(myData);

            if (iter.good())
            {
                label nbrBFacei = *iter;

                masterCoupledFaces[coupledI] = mesh.nInternalFaces() + bFacei;
                slaveCoupledFaces[coupledI] =
                    domainMesh.nInternalFaces()
                  + nbrBFacei;

                coupledI++;
            }
        }
    }

    masterCoupledFaces.setSize(coupledI);
    slaveCoupledFaces.setSize(coupledI);

    if (debug)
    {
        Pout<< "findCouples : found " << coupledI
            << " faces that will be stitched" << nl << endl;
    }
}


void Foam::fvMeshDistribute::findCouples
(
    const UPtrList<polyMesh>& meshes,
    const PtrList<labelList>& domainSourceFaces,
    const PtrList<labelList>& domainSourceProcs,
    const PtrList<labelList>& domainSourcePatchs,

    labelListList& localBoundaryFace,
    labelListList& remoteFaceProc,
    labelListList& remoteBoundaryFace
)
{
    // Collect all matching processor face pairs. These are all the
    // faces for which we have both sides

    // Pass 0: count number of inter-processor faces
    labelList nProcFaces(meshes.size(), 0);
    forAll(meshes, meshi)
    {
        if (meshes.set(meshi))
        {
            const labelList& domainProc = domainSourceProcs[meshi];
            const labelList& domainPatch = domainSourcePatchs[meshi];

            forAll(domainProc, bFacei)
            {
                if (domainProc[bFacei] != -1 && domainPatch[bFacei] == -1)
                {
                    nProcFaces[meshi]++;
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "fvMeshDistribute::findCouples : nProcFaces:"
            << flatOutput(nProcFaces) << endl;
    }


    // Size
    List<DynamicList<label>> dynLocalFace(Pstream::nProcs());
    List<DynamicList<label>> dynRemoteProc(Pstream::nProcs());
    List<DynamicList<label>> dynRemoteFace(Pstream::nProcs());

    forAll(meshes, meshi)
    {
        if (meshes.set(meshi))
        {
            dynLocalFace[meshi].reserve(nProcFaces[meshi]);
            dynRemoteProc[meshi].reserve(nProcFaces[meshi]);
            dynRemoteFace[meshi].reserve(nProcFaces[meshi]);
        }
    }


    // Insert all into big map. Find matches
    LabelPairMap<labelPair> map(2*sum(nProcFaces));

    nProcFaces = 0;

    forAll(meshes, meshi)
    {
        if (meshes.set(meshi))
        {
            const labelList& domainFace = domainSourceFaces[meshi];
            const labelList& domainProc = domainSourceProcs[meshi];
            const labelList& domainPatch = domainSourcePatchs[meshi];

            forAll(domainProc, bFacei)
            {
                if (domainProc[bFacei] != -1 && domainPatch[bFacei] == -1)
                {
                    const labelPair key
                    (
                        domainProc[bFacei],
                        domainFace[bFacei]
                    );
                    auto fnd = map.find(key);

                    if (!fnd.good())
                    {
                        // Insert
                        map.emplace(key, meshi, bFacei);
                    }
                    else
                    {
                        // Second occurence. Found match.
                        const label matchProci = fnd().first();
                        const label matchFacei = fnd().second();

                        dynLocalFace[meshi].append(bFacei);
                        dynRemoteProc[meshi].append(matchProci);
                        dynRemoteFace[meshi].append(matchFacei);
                        nProcFaces[meshi]++;

                        dynLocalFace[matchProci].append(matchFacei);
                        dynRemoteProc[matchProci].append(meshi);
                        dynRemoteFace[matchProci].append(bFacei);
                        nProcFaces[matchProci]++;
                    }
                }
            }
        }
    }

    if (debug)
    {
        Pout<< "fvMeshDistribute::findCouples : stored procFaces:"
            << map.size() << endl;
    }

    localBoundaryFace.setSize(Pstream::nProcs());
    remoteFaceProc.setSize(Pstream::nProcs());
    remoteBoundaryFace.setSize(Pstream::nProcs());
    forAll(meshes, meshi)
    {
        if (meshes.set(meshi))
        {
            localBoundaryFace[meshi] = std::move(dynLocalFace[meshi]);
            remoteFaceProc[meshi] = std::move(dynRemoteProc[meshi]);
            remoteBoundaryFace[meshi] = std::move(dynRemoteFace[meshi]);
        }
    }


    if (debug)
    {
        Pout<< "fvMeshDistribute::findCouples : found matches:"
            << flatOutput(nProcFaces) << endl;
    }
}


// Map data on boundary faces to new mesh (resulting from adding two meshes)
Foam::labelList Foam::fvMeshDistribute::mapBoundaryData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const label nInternalFaces1,
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nBoundaryFaces());

    forAll(boundaryData0, oldBFacei)
    {
        label newFacei = map.oldFaceMap()[oldBFacei + map.nOldInternalFaces()];

        // Face still exists (is necessary?) and still boundary face
        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData0[oldBFacei];
        }
    }

    forAll(boundaryData1, addedBFacei)
    {
        label newFacei = map.addedFaceMap()[addedBFacei + nInternalFaces1];

        if (newFacei >= 0 && newFacei >= mesh.nInternalFaces())
        {
            newBoundaryData[newFacei - mesh.nInternalFaces()] =
                boundaryData1[addedBFacei];
        }
    }

    return newBoundaryData;
}


Foam::labelList Foam::fvMeshDistribute::mapPointData
(
    const primitiveMesh& mesh,      // mesh after adding
    const mapAddedPolyMesh& map,
    const labelList& boundaryData0, // on mesh before adding
    const labelList& boundaryData1  // on added mesh
)
{
    labelList newBoundaryData(mesh.nPoints());

    forAll(boundaryData0, oldPointi)
    {
        label newPointi = map.oldPointMap()[oldPointi];

        // Point still exists (is necessary?)
        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData0[oldPointi];
        }
    }

    forAll(boundaryData1, addedPointi)
    {
        label newPointi = map.addedPointMap()[addedPointi];

        if (newPointi >= 0)
        {
            newBoundaryData[newPointi] = boundaryData1[addedPointi];
        }
    }

    return newBoundaryData;
}


// Remove cells. Add all exposed faces to patch oldInternalPatchi
Foam::autoPtr<Foam::mapPolyMesh> Foam::fvMeshDistribute::doRemoveCells
(
    const labelList& cellsToRemove,
    const label oldInternalPatchi
)
{
    // Mesh change engine
    polyTopoChange meshMod(mesh_);

    // Cell removal topo engine. Do NOT synchronize parallel since
    // we are doing a local cell removal.
    removeCells cellRemover(mesh_, false);

    // Get all exposed faces
    labelList exposedFaces(cellRemover.getExposedFaces(cellsToRemove));

    // Insert the topo changes
    cellRemover.setRefinement
    (
        cellsToRemove,
        exposedFaces,
        labelList(exposedFaces.size(), oldInternalPatchi),  // patch for exposed
                                                            // faces.
        meshMod
    );


    //// Generate test field
    //tmp<surfaceScalarField> sfld(generateTestField(mesh_));

    // Save internal fields (note: not as DimensionedFields since would
    // get mapped)
    PtrList<Field<scalar>> sFlds;
    saveInternalFields(sFlds);
    PtrList<Field<vector>> vFlds;
    saveInternalFields(vFlds);
    PtrList<Field<sphericalTensor>> sptFlds;
    saveInternalFields(sptFlds);
    PtrList<Field<symmTensor>> sytFlds;
    saveInternalFields(sytFlds);
    PtrList<Field<tensor>> tFlds;
    saveInternalFields(tFlds);

    // Change the mesh. No inflation. Note: no parallel comms allowed.

    // TBD: temporarily unset mesh moving to avoid problems in meshflux
    //      mapping. To be investigated.
    const bool oldMoving = mesh_.moving(false);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(mesh_, false, false);
    mesh_.moving(oldMoving);

    // Update fields
    mesh_.updateMesh(map());


    // Any exposed faces in a surfaceField will not be mapped. Map the value
    // of these separately (until there is support in all PatchFields for
    // mapping from internal faces ...)

    mapExposedFaces(map(), sFlds);
    mapExposedFaces(map(), vFlds);
    mapExposedFaces(map(), sptFlds);
    mapExposedFaces(map(), sytFlds);
    mapExposedFaces(map(), tFlds);


    //// Test test field
    //testField(sfld);


    // Move mesh (since morphing does not do this)
    if (map().hasMotionPoints())
    {
        mesh_.movePoints(map().preMotionPoints());
    }


    return map;
}


// Delete and add processor patches. Changes mesh and returns per neighbour proc
// the processor patchID.
void Foam::fvMeshDistribute::addProcPatches
(
    const labelList& nbrProc,       // processor that neighbour is now on
    const labelList& referPatchID,  // patchID (or -1) I originated from
    List<Map<label>>& procPatchID
)
{
    // Now use the neighbourFace/Proc to repatch the mesh. These lists
    // contain for all current boundary faces the global patchID (for non-proc
    // patch) or the processor.

    // Determine a visit order such that the processor patches get added
    // in order of increasing neighbour processor (and for same neighbour
    // processor (in case of processor cyclics) in order of increasing
    // 'refer' patch)
    labelList indices;
    sortedOrder(nbrProc, indices, lessProcPatches(nbrProc, referPatchID));

    procPatchID.setSize(Pstream::nProcs());

    forAll(indices, i)
    {
        label bFacei = indices[i];
        label proci = nbrProc[bFacei];

        if (proci != -1 && proci != Pstream::myProcNo())
        {
            if (!procPatchID[proci].found(referPatchID[bFacei]))
            {
                // No patch for neighbour yet. Is either a normal processor
                // patch or a processorCyclic patch.

                if (referPatchID[bFacei] == -1)
                {
                    // Ordinary processor boundary

                    processorPolyPatch pp
                    (
                        0,              // size
                        mesh_.nFaces(),
                        mesh_.boundaryMesh().size(),
                        mesh_.boundaryMesh(),
                        Pstream::myProcNo(),
                        proci
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh_,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
                else
                {
                    const coupledPolyPatch& pcPatch
                        = refCast<const coupledPolyPatch>
                          (
                              mesh_.boundaryMesh()[referPatchID[bFacei]]
                          );
                    processorCyclicPolyPatch pp
                    (
                        0,              // size
                        mesh_.nFaces(),
                        mesh_.boundaryMesh().size(),
                        mesh_.boundaryMesh(),
                        Pstream::myProcNo(),
                        proci,
                        pcPatch.name(),
                        pcPatch.transform()
                    );

                    procPatchID[proci].insert
                    (
                        referPatchID[bFacei],
                        fvMeshTools::addPatch
                        (
                            mesh_,
                            pp,
                            dictionary(),   // optional per field patchField
                            processorCyclicFvPatchField<scalar>::typeName,
                            false           // not parallel sync
                        )
                    );
                }
            }
        }
    }
}


// Get boundary faces to be repatched. Is -1 or new patchID
Foam::labelList Foam::fvMeshDistribute::getBoundaryPatch
(
    const labelList& nbrProc,               // new processor per boundary face
    const labelList& referPatchID,          // patchID (or -1) I originated from
    const List<Map<label>>& procPatchID    // per proc the new procPatches
)
{
    labelList patchIDs(nbrProc);

    forAll(nbrProc, bFacei)
    {
        if (nbrProc[bFacei] == Pstream::myProcNo())
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = origPatchi;
        }
        else if (nbrProc[bFacei] != -1)
        {
            label origPatchi = referPatchID[bFacei];
            patchIDs[bFacei] = procPatchID[nbrProc[bFacei]][origPatchi];
        }
        else
        {
            patchIDs[bFacei] = -1;
        }
    }
    return patchIDs;
}


// Send mesh and coupling data.
void Foam::fvMeshDistribute::sendMesh
(
    const label domain,
    const fvMesh& mesh,

    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,

    const labelList& sourceFace,
    const labelList& sourceProc,
    const labelList& sourcePatch,
    const labelList& sourceNewNbrProc,
    const labelList& sourcePointMaster,
    Ostream& toDomain
)
{
    if (debug)
    {
        Pout<< "Sending to domain " << domain << nl
            << "    nPoints:" << mesh.nPoints() << nl
            << "    nFaces:" << mesh.nFaces() << nl
            << "    nCells:" << mesh.nCells() << nl
            << "    nPatches:" << mesh.boundaryMesh().size() << nl
            << endl;
    }

    // Assume sparse, possibly overlapping point zones. Get contents
    // in merged-zone indices.
    CompactListList<label> zonePoints;
    {
        const pointZoneMesh& pointZones = mesh.pointZones();

        labelList rowSizes(pointZoneNames.size(), Zero);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = pointZones[myZoneID].size();
            }
        }
        zonePoints.setSize(rowSizes);

        forAll(pointZoneNames, nameI)
        {
            label myZoneID = pointZones.findZoneID(pointZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zonePoints[nameI].deepCopy(pointZones[myZoneID]);
            }
        }
    }

    // Assume sparse, possibly overlapping face zones
    CompactListList<label> zoneFaces;
    CompactListList<bool> zoneFaceFlip;
    {
        const faceZoneMesh& faceZones = mesh.faceZones();

        labelList rowSizes(faceZoneNames.size(), Zero);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = faceZones[myZoneID].size();
            }
        }

        zoneFaces.setSize(rowSizes);
        zoneFaceFlip.setSize(rowSizes);

        forAll(faceZoneNames, nameI)
        {
            label myZoneID = faceZones.findZoneID(faceZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneFaces[nameI].deepCopy(faceZones[myZoneID]);
                zoneFaceFlip[nameI].deepCopy(faceZones[myZoneID].flipMap());
            }
        }
    }

    // Assume sparse, possibly overlapping cell zones
    CompactListList<label> zoneCells;
    {
        const cellZoneMesh& cellZones = mesh.cellZones();

        labelList rowSizes(cellZoneNames.size(), Zero);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                rowSizes[nameI] = cellZones[myZoneID].size();
            }
        }

        zoneCells.setSize(rowSizes);

        forAll(cellZoneNames, nameI)
        {
            label myZoneID = cellZones.findZoneID(cellZoneNames[nameI]);

            if (myZoneID != -1)
            {
                zoneCells[nameI].deepCopy(cellZones[myZoneID]);
            }
        }
    }
    ////- Assume full cell zones
    //labelList cellZoneID;
    //if (hasCellZones)
    //{
    //    cellZoneID.setSize(mesh.nCells());
    //    cellZoneID = -1;
    //
    //    const cellZoneMesh& cellZones = mesh.cellZones();
    //
    //    forAll(cellZones, zoneI)
    //    {
    //        labelUIndList(cellZoneID, cellZones[zoneI]) = zoneI;
    //    }
    //}

    // Send
    toDomain
        << mesh.points()
        // Send as (offsets/values)
        << CompactListList<label>::pack<face>(mesh.faces())
        << mesh.faceOwner()
        << mesh.faceNeighbour()
        << mesh.boundaryMesh()

        << zonePoints
        << zoneFaces
        << zoneFaceFlip
        << zoneCells

        << sourceFace
        << sourceProc
        << sourcePatch
        << sourceNewNbrProc
        << sourcePointMaster;


    if (debug)
    {
        Pout<< "Started sending mesh to domain " << domain
            << endl;
    }
}


// Receive mesh. Opposite of sendMesh
Foam::autoPtr<Foam::fvMesh> Foam::fvMeshDistribute::receiveMesh
(
    const label domain,
    const wordList& pointZoneNames,
    const wordList& faceZoneNames,
    const wordList& cellZoneNames,
    const Time& runTime,
    labelList& domainSourceFace,
    labelList& domainSourceProc,
    labelList& domainSourcePatch,
    labelList& domainSourceNewNbrProc,
    labelList& domainSourcePointMaster,
    Istream& fromNbr
)
{
    pointField domainPoints(fromNbr);
    // Receive as (offsets/values), unpack to faceList
    faceList domainFaces = CompactListList<label>(fromNbr).unpack<face>();
    labelList domainAllOwner(fromNbr);
    labelList domainAllNeighbour(fromNbr);
    PtrList<entry> patchEntries(fromNbr);

    CompactListList<label> zonePoints(fromNbr);
    CompactListList<label> zoneFaces(fromNbr);
    CompactListList<bool> zoneFaceFlip(fromNbr);
    CompactListList<label> zoneCells(fromNbr);

    fromNbr
        >> domainSourceFace
        >> domainSourceProc
        >> domainSourcePatch
        >> domainSourceNewNbrProc
        >> domainSourcePointMaster;

    // Construct fvMesh
    auto domainMeshPtr = autoPtr<fvMesh>::New
    (
        IOobject
        (
            polyMesh::defaultRegion,
            runTime.timeName(),
            runTime,
            IOobject::NO_READ
        ),
        std::move(domainPoints),
        std::move(domainFaces),
        std::move(domainAllOwner),
        std::move(domainAllNeighbour),
        false                   // no parallel comms
    );
    fvMesh& domainMesh = *domainMeshPtr;

    polyPatchList patches(patchEntries.size());

    forAll(patchEntries, patchi)
    {
        patches.set
        (
            patchi,
            polyPatch::New
            (
                patchEntries[patchi].keyword(),
                patchEntries[patchi].dict(),
                patchi,
                domainMesh.boundaryMesh()
            )
        );
    }
    // Add patches; no parallel comms
    domainMesh.addFvPatches(patches, false);

    // Construct zones
    List<pointZone*> pZonePtrs(pointZoneNames.size());
    forAll(pZonePtrs, i)
    {
        pZonePtrs[i] = new pointZone
        (
            pointZoneNames[i],
            zonePoints[i],
            i,
            domainMesh.pointZones()
        );
    }

    List<faceZone*> fZonePtrs(faceZoneNames.size());
    forAll(fZonePtrs, i)
    {
        fZonePtrs[i] = new faceZone
        (
            faceZoneNames[i],
            zoneFaces[i],
            zoneFaceFlip[i],
            i,
            domainMesh.faceZones()
        );
    }

    List<cellZone*> cZonePtrs(cellZoneNames.size());
    forAll(cZonePtrs, i)
    {
        cZonePtrs[i] = new cellZone
        (
            cellZoneNames[i],
            zoneCells[i],
            i,
            domainMesh.cellZones()
        );
    }
    domainMesh.addZones(pZonePtrs, fZonePtrs, cZonePtrs);

    return domainMeshPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fvMeshDistribute::fvMeshDistribute(fvMesh& mesh)
:
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::fvMeshDistribute::countCells
(
    const labelList& distribution
)
{
    labelList nCells(Pstream::nProcs(), Zero);
    forAll(distribution, celli)
    {
        label newProc = distribution[celli];

        if (newProc < 0 || newProc >= Pstream::nProcs())
        {
            FatalErrorInFunction
                << "Distribution should be in range 0.." << Pstream::nProcs()-1
                << endl
                << "At index " << celli << " distribution:" << newProc
                << abort(FatalError);
        }
        nCells[newProc]++;
    }
    return nCells;
}


Foam::autoPtr<Foam::mapDistributePolyMesh> Foam::fvMeshDistribute::distribute
(
    const labelList& distribution
)
{
    // Some checks on distribution
    if (distribution.size() != mesh_.nCells())
    {
        FatalErrorInFunction
            << "Size of distribution:"
            << distribution.size() << " mesh nCells:" << mesh_.nCells()
            << abort(FatalError);
    }


    const polyBoundaryMesh& patches = mesh_.boundaryMesh();

    // Check all processors have same non-proc patches in same order.
    if (patches.checkParallelSync(true))
    {
        FatalErrorInFunction
            << "This application requires all non-processor patches"
            << " to be present in the same order on all patches" << nl
            << "followed by the processor patches (which of course are unique)."
            << nl
            << "Local patches:" << mesh_.boundaryMesh().names()
            << abort(FatalError);
    }

    // Save some data for mapping later on
    const label nOldPoints(mesh_.nPoints());
    const label nOldFaces(mesh_.nFaces());
    const label nOldCells(mesh_.nCells());
    labelList oldPatchStarts(patches.size());
    labelList oldPatchNMeshPoints(patches.size());
    forAll(patches, patchi)
    {
        oldPatchStarts[patchi] = patches[patchi].start();
        oldPatchNMeshPoints[patchi] = patches[patchi].nPoints();
    }



    // Short circuit trivial case.
    if (!Pstream::parRun())
    {
        // Collect all maps and return
        return autoPtr<mapDistributePolyMesh>::New
        (
            mesh_,

            nOldPoints,
            nOldFaces,
            nOldCells,
            std::move(oldPatchStarts),
            std::move(oldPatchNMeshPoints),

            labelListList(one(), identity(mesh_.nPoints())), //subPointMap
            labelListList(one(), identity(mesh_.nFaces())),  //subFaceMap
            labelListList(one(), identity(mesh_.nCells())),  //subCellMap
            labelListList(one(), identity(patches.size())),  //subPatchMap

            labelListList(one(), identity(mesh_.nPoints())), //pointMap
            labelListList(one(), identity(mesh_.nFaces())),  //faceMap
            labelListList(one(), identity(mesh_.nCells())),  //cellMap
            labelListList(one(), identity(patches.size()))   //patchMap
        );
    }


    // Collect any zone names over all processors and shuffle zones accordingly
    // Note that this is not necessary for redistributePar since that already
    // checks for it. However other use (e.g. mesh generators) might not.
    const wordList pointZoneNames(mergeWordList(mesh_.pointZones().names()));
    reorderZones<pointZone>(pointZoneNames, mesh_.pointZones());

    const wordList faceZoneNames(mergeWordList(mesh_.faceZones().names()));
    reorderZones<faceZone>(faceZoneNames, mesh_.faceZones());

    const wordList cellZoneNames(mergeWordList(mesh_.cellZones().names()));
    reorderZones<cellZone>(cellZoneNames, mesh_.cellZones());


    // Local environment of all boundary faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // A face is uniquely defined by
    //  - proc
    //  - local face no
    //
    // To glue the parts of meshes which can get sent from anywhere we
    // need to know on boundary faces what the above tuple on both sides is.
    // So we need to maintain:
    //  - original face
    //  - original processor id (= trivial)
    // For coupled boundaries (where the faces are 'duplicate') we take the
    // lowest numbered processor as the data to store.
    //
    // Additionally to create the procboundaries we need to know where the owner
    // cell on the other side now is: newNeighbourProc.
    //

    // physical boundary:
    //     sourceProc = -1
    //     sourceNewNbrProc = -1
    //     sourceFace = -1
    //     sourcePatch = patchID
    // processor boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = -1
    // ?cyclic:
    // ?    sourceProc = proc
    // ?    sourceNewNbrProc = distribution of coupled cell
    // ?    sourceFace = face (on owner side)
    // ?    sourcePatch = patchID
    // processor-cyclic boundary:
    //     sourceProc = proc (on owner side)
    //     sourceNewNbrProc = distribution of coupled cell
    //     sourceFace = face (on owner side)
    //     sourcePatch = patchID

    labelList sourcePatch;
    labelList sourceFace;
    labelList sourceProc;
    labelList sourceNewNbrProc;
    labelList sourcePointMaster;
    getCouplingData
    (
        distribution,
        sourceFace,
        sourceProc,
        sourcePatch,
        sourceNewNbrProc,
        sourcePointMaster
    );


    // Remove meshPhi. Since this would otherwise disappear anyway
    // during topo changes and we have to guarantee that all the fields
    // can be sent.

    // NOTE: could/should use (isMeshUpdate = true) for mesh_.clearOut()
    // but the bottom level will do a clearGeom() and that doesn't seem
    // to work particularly well with isMeshUpdate at all.
    // re-visit if needed...

    mesh_.clearOut();

    mesh_.resetMotion();

    // Get data to send. Make sure is synchronised

    HashTable<wordList> allFieldNames;

    getFieldNames<volScalarField>(mesh_, allFieldNames);
    getFieldNames<volVectorField>(mesh_, allFieldNames);
    getFieldNames<volSphericalTensorField>(mesh_, allFieldNames);
    getFieldNames<volSymmTensorField>(mesh_, allFieldNames);
    getFieldNames<volTensorField>(mesh_, allFieldNames);

    getFieldNames<surfaceScalarField>(mesh_, allFieldNames);
    getFieldNames<surfaceVectorField>(mesh_, allFieldNames);
    getFieldNames<surfaceSphericalTensorField>(mesh_, allFieldNames);
    getFieldNames<surfaceSymmTensorField>(mesh_, allFieldNames);
    getFieldNames<surfaceTensorField>(mesh_, allFieldNames);

    getFieldNames<volScalarField::Internal>
    (
        mesh_,
        allFieldNames,
        volScalarField::typeName
    );
    getFieldNames<volVectorField::Internal>
    (
        mesh_,
        allFieldNames,
        volVectorField::typeName
    );
    getFieldNames<volSphericalTensorField::Internal>
    (
        mesh_,
        allFieldNames,
        volSphericalTensorField::typeName
    );
    getFieldNames<volSymmTensorField::Internal>
    (
        mesh_,
        allFieldNames,
        volSymmTensorField::typeName
    );
    getFieldNames<volTensorField::Internal>
    (
        mesh_,
        allFieldNames,
        volTensorField::typeName
    );


    // Find patch to temporarily put exposed and processor faces into.
    const label oldInternalPatchi = findNonEmptyPatch();


    // Delete processor patches, starting from the back. Move all faces into
    // oldInternalPatchi.
    labelList repatchFaceMap;
    {
        autoPtr<mapPolyMesh> repatchMap = deleteProcPatches(oldInternalPatchi);

        // Store face map (only face ordering that changed)
        repatchFaceMap = repatchMap().faceMap();

        // Reorder all boundary face data (sourceProc, sourceFace etc.)
        labelList bFaceMap
        (
            SubList<label>
            (
                repatchMap().reverseFaceMap(),
                mesh_.nBoundaryFaces(),
                mesh_.nInternalFaces()
            )
          - mesh_.nInternalFaces()
        );

        inplaceReorder(bFaceMap, sourceFace);
        inplaceReorder(bFaceMap, sourceProc);
        inplaceReorder(bFaceMap, sourcePatch);
        inplaceReorder(bFaceMap, sourceNewNbrProc);
    }



    // Print a bit.
    if (debug)
    {
        Pout<< nl << "MESH WITH PROC PATCHES DELETED:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printIntFieldInfo<volScalarField::Internal>(mesh_);
        printIntFieldInfo<volVectorField::Internal>(mesh_);
        printIntFieldInfo<volSphericalTensorField::Internal>(mesh_);
        printIntFieldInfo<volSymmTensorField::Internal>(mesh_);
        printIntFieldInfo<volTensorField::Internal>(mesh_);
        Pout<< nl << endl;
    }



    // Maps from subsetted mesh (that is sent) back to original maps
    labelListList subCellMap(Pstream::nProcs());
    labelListList subFaceMap(Pstream::nProcs());
    labelListList subPointMap(Pstream::nProcs());
    labelListList subPatchMap(Pstream::nProcs());
    // Maps from subsetted mesh to reconstructed mesh
    labelListList constructCellMap(Pstream::nProcs());
    labelListList constructFaceMap(Pstream::nProcs());
    labelListList constructPointMap(Pstream::nProcs());
    labelListList constructPatchMap(Pstream::nProcs());


    // Find out schedule
    // ~~~~~~~~~~~~~~~~~

    labelList nSendCells(countCells(distribution));
    labelList nRecvCells(Pstream::nProcs());
    UPstream::allToAll(nSendCells, nRecvCells);

    // Allocate buffers
    PstreamBuffers pBufs;


    // What to send to neighbouring domains
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Disable parallel.
    const bool oldParRun = UPstream::parRun(false);

    forAll(nSendCells, recvProc)
    {
        if (recvProc != Pstream::myProcNo() && nSendCells[recvProc] > 0)
        {
            // Send to recvProc

            if (debug)
            {
                Pout<< nl
                    << "SUBSETTING FOR DOMAIN " << recvProc
                    << " cells to send:"
                    << nSendCells[recvProc]
                    << nl << endl;
            }

            // Pstream for sending mesh and fields
            UOPstream str(recvProc, pBufs);

            // Mesh subsetting engine - subset the cells of the current domain.
            fvMeshSubset subsetter
            (
                mesh_,
                recvProc,
                distribution,
                oldInternalPatchi,  // oldInternalFaces patch
                false               // no parallel sync
            );

            subCellMap[recvProc] = subsetter.cellMap();
            subFaceMap[recvProc] = subsetter.faceFlipMap();
            inplaceRenumberWithFlip
            (
                repatchFaceMap,
                false,      // oldToNew has flip
                true,       // subFaceMap has flip
                subFaceMap[recvProc]
            );
            subPointMap[recvProc] = subsetter.pointMap();
            subPatchMap[recvProc] = subsetter.patchMap();


            // Subset the boundary fields (owner/neighbour/processor)
            labelList procSourceFace;
            labelList procSourceProc;
            labelList procSourcePatch;
            labelList procSourceNewNbrProc;
            labelList procSourcePointMaster;

            subsetCouplingData
            (
                subsetter.subMesh(),
                subsetter.pointMap(),       // from subMesh to mesh
                subsetter.faceMap(),        //      ,,      ,,
                subsetter.cellMap(),        //      ,,      ,,

                distribution,               // old mesh distribution
                mesh_.faceOwner(),          // old owner
                mesh_.faceNeighbour(),
                mesh_.nInternalFaces(),

                sourceFace,
                sourceProc,
                sourcePatch,
                sourceNewNbrProc,
                sourcePointMaster,

                procSourceFace,
                procSourceProc,
                procSourcePatch,
                procSourceNewNbrProc,
                procSourcePointMaster
            );


            // Send to neighbour
            sendMesh
            (
                recvProc,
                subsetter.subMesh(),

                pointZoneNames,
                faceZoneNames,
                cellZoneNames,

                procSourceFace,
                procSourceProc,
                procSourcePatch,
                procSourceNewNbrProc,
                procSourcePointMaster,

                str
            );

            // volFields
            sendFields<volScalarField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volVectorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volSphericalTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volSymmTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );

            // surfaceFields
            sendFields<surfaceScalarField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<surfaceVectorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<surfaceSphericalTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<surfaceSymmTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<surfaceTensorField>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );

            // Dimensioned fields
            sendFields<volScalarField::Internal>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volVectorField::Internal>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volSphericalTensorField::Internal>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volSymmTensorField::Internal>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
            sendFields<volTensorField::Internal>
            (
                recvProc,
                allFieldNames,
                subsetter,
                str
            );
        }
    }


    UPstream::parRun(oldParRun);  // Restore parallel state

    if (debug)
    {
        Pout<< "Starting sending" << endl;
    }

    pBufs.finishedSends();

    if (debug)
    {
        Pout<< "Finished sending and receiving : "
            << flatOutput(pBufs.recvDataCounts()) << endl;
    }


    // Subset the part that stays
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        // Save old mesh maps before changing mesh
        const labelList oldFaceOwner(mesh_.faceOwner());
        const labelList oldFaceNeighbour(mesh_.faceNeighbour());
        const label oldInternalFaces = mesh_.nInternalFaces();

        // Remove cells.
        autoPtr<mapPolyMesh> subMap
        (
            doRemoveCells
            (
                select(false, distribution, Pstream::myProcNo()),
                oldInternalPatchi
            )
        );

        // Addressing from subsetted mesh
        subCellMap[Pstream::myProcNo()] = subMap().cellMap();
        subFaceMap[Pstream::myProcNo()] = renumber
        (
            repatchFaceMap,
            subMap().faceMap()
        );
        // Insert the sign bit from face flipping
        labelList& faceMap = subFaceMap[Pstream::myProcNo()];
        forAll(faceMap, faceI)
        {
            faceMap[faceI] += 1;
        }
        const labelHashSet& flip = subMap().flipFaceFlux();
        for (const label facei : flip)
        {
            faceMap[facei] = -faceMap[facei];
        }
        subPointMap[Pstream::myProcNo()] = subMap().pointMap();
        subPatchMap[Pstream::myProcNo()] = identity(patches.size());

        // Subset the mesh data: neighbourCell/neighbourProc fields
        labelList domainSourceFace;
        labelList domainSourceProc;
        labelList domainSourcePatch;
        labelList domainSourceNewNbrProc;
        labelList domainSourcePointMaster;

        subsetCouplingData
        (
            mesh_,                          // new mesh
            subMap().pointMap(),            // from new to original mesh
            subMap().faceMap(),             // from new to original mesh
            subMap().cellMap(),

            distribution,                   // distribution before subsetting
            oldFaceOwner,                   // owner before subsetting
            oldFaceNeighbour,               // neighbour        ,,
            oldInternalFaces,               // nInternalFaces   ,,

            sourceFace,
            sourceProc,
            sourcePatch,
            sourceNewNbrProc,
            sourcePointMaster,

            domainSourceFace,
            domainSourceProc,
            domainSourcePatch,
            domainSourceNewNbrProc,
            domainSourcePointMaster
        );

        sourceFace.transfer(domainSourceFace);
        sourceProc.transfer(domainSourceProc);
        sourcePatch.transfer(domainSourcePatch);
        sourceNewNbrProc.transfer(domainSourceNewNbrProc);
        sourcePointMaster.transfer(domainSourcePointMaster);
    }


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "STARTING MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printIntFieldInfo<volScalarField::Internal>(mesh_);
        printIntFieldInfo<volVectorField::Internal>(mesh_);
        printIntFieldInfo<volSphericalTensorField::Internal>(mesh_);
        printIntFieldInfo<volSymmTensorField::Internal>(mesh_);
        printIntFieldInfo<volTensorField::Internal>(mesh_);
        Pout<< nl << endl;
    }



    // Receive and add what was sent
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Disable parallel. Original state already known.
    UPstream::parRun(false);

    PtrList<labelList> domainSourceFaces(Pstream::nProcs());
    PtrList<labelList> domainSourceProcs(Pstream::nProcs());
    PtrList<labelList> domainSourcePatchs(Pstream::nProcs());
    PtrList<labelList> domainSourceNewNbrProcs(Pstream::nProcs());
    PtrList<labelList> domainSourcePointMasters(Pstream::nProcs());

    PtrList<fvMesh> domainMeshPtrs(Pstream::nProcs());

    // Define field storage lists
    #undef  doLocalCode
    #define doLocalCode(FieldType, Variable)                                  \
        PtrList<PtrList<FieldType>> Variable(UPstream::nProcs());

    doLocalCode(volScalarField, vsfs);
    doLocalCode(volVectorField, vvfs);
    doLocalCode(volSphericalTensorField, vsptfs);
    doLocalCode(volSymmTensorField, vsytfs);
    doLocalCode(volTensorField, vtfs);

    doLocalCode(surfaceScalarField, ssfs);
    doLocalCode(surfaceVectorField, svfs);
    doLocalCode(surfaceSphericalTensorField, ssptfs);
    doLocalCode(surfaceSymmTensorField, ssytfs);
    doLocalCode(surfaceTensorField, stfs);

    doLocalCode(volScalarField::Internal, dsfs);
    doLocalCode(volVectorField::Internal, dvfs);
    doLocalCode(volSphericalTensorField::Internal, dsptfs)
    doLocalCode(volSymmTensorField::Internal, dsytfs);
    doLocalCode(volTensorField::Internal, dtfs);

    #undef doLocalCode


    forAll(nRecvCells, sendProc)
    {
        // Did processor sendProc send anything to me?
        if (sendProc != UPstream::myProcNo() && nRecvCells[sendProc] > 0)
        {
            if (debug)
            {
                Pout<< nl
                    << "RECEIVING FROM DOMAIN " << sendProc
                    << " cells to receive:"
                    << nRecvCells[sendProc]
                    << nl << endl;
            }


            // Pstream for receiving mesh and fields
            UIPstream str(sendProc, pBufs);


            // Receive from sendProc - opposite of sendMesh
            {
                autoPtr<fvMesh> domainMeshPtr = receiveMesh
                (
                    sendProc,
                    pointZoneNames,
                    faceZoneNames,
                    cellZoneNames,

                    const_cast<Time&>(mesh_.time()),

                    domainSourceFaces.emplace_set(sendProc),
                    domainSourceProcs.emplace_set(sendProc),
                    domainSourcePatchs.emplace_set(sendProc),
                    domainSourceNewNbrProcs.emplace_set(sendProc),
                    domainSourcePointMasters.emplace_set(sendProc),
                    str
                );
                domainMeshPtrs.set(sendProc, std::move(domainMeshPtr));
                fvMesh& domainMesh = domainMeshPtrs[sendProc];
                // Force construction of various on mesh.
                //(void)domainMesh.globalData();


                // Receive fields. Read as single dictionary because
                // of problems reading consecutive fields from single stream.
                dictionary fieldDicts(str);

                #undef  doLocalCode
                #define doLocalCode(FieldType, Variable)                      \
                    receiveFields<FieldType>                                  \
                    (                                                         \
                        sendProc,                                             \
                        allFieldNames,                                        \
                        domainMesh,                                           \
                        Variable.emplace_set(sendProc),                       \
                        fieldDicts                                            \
                    )

                // Volume Fields
                doLocalCode(volScalarField, vsfs);
                doLocalCode(volVectorField, vvfs);
                doLocalCode(volSphericalTensorField, vsptfs);
                doLocalCode(volSymmTensorField, vsytfs);
                doLocalCode(volTensorField, vtfs);

                // Surface Fields
                doLocalCode(surfaceScalarField, ssfs);
                doLocalCode(surfaceVectorField, svfs);
                doLocalCode(surfaceSphericalTensorField, ssptfs);
                doLocalCode(surfaceSymmTensorField, ssytfs);
                doLocalCode(surfaceTensorField, stfs);

                // Dimensioned Fields
                doLocalCode(volScalarField::Internal, dsfs);
                doLocalCode(volVectorField::Internal, dvfs);
                doLocalCode(volSphericalTensorField::Internal, dsptfs);
                doLocalCode(volSymmTensorField::Internal, dsytfs);
                doLocalCode(volTensorField::Internal, dtfs);

                #undef doLocalCode
            }
        }
    }

    // Clear out storage
    pBufs.clear();


    // Set up pointers to meshes so we can include our mesh_
    UPtrList<polyMesh> meshes(domainMeshPtrs.size());
    UPtrList<fvMesh> fvMeshes(domainMeshPtrs.size());
    forAll(domainMeshPtrs, proci)
    {
        if (domainMeshPtrs.set(proci))
        {
            meshes.set(proci, &domainMeshPtrs[proci]);
            fvMeshes.set(proci, &domainMeshPtrs[proci]);
        }
    }

    // 'Receive' from myself
    {
        meshes.set(UPstream::myProcNo(), &mesh_);
        fvMeshes.set(UPstream::myProcNo(), &mesh_);

        domainSourceFaces.emplace_set(UPstream::myProcNo()) = sourceFace;
        domainSourceProcs.emplace_set(UPstream::myProcNo()) = sourceProc;
        domainSourcePatchs.emplace_set(UPstream::myProcNo()) = sourcePatch;

        domainSourceNewNbrProcs.emplace_set(UPstream::myProcNo()) =
            sourceNewNbrProc;

        domainSourcePointMasters.emplace_set(UPstream::myProcNo()) =
            sourcePointMaster;
    }


    // Find matching faces that need to be stitched
    labelListList localBoundaryFace(UPstream::nProcs());
    labelListList remoteFaceProc(UPstream::nProcs());
    labelListList remoteBoundaryFace(UPstream::nProcs());
    findCouples
    (
        meshes,
        domainSourceFaces,
        domainSourceProcs,
        domainSourcePatchs,

        localBoundaryFace,
        remoteFaceProc,
        remoteBoundaryFace
    );


    const label nOldInternalFaces = mesh_.nInternalFaces();
    const labelList oldFaceOwner(mesh_.faceOwner());

    // TBD: temporarily unset mesh moving to avoid problems in meshflux
    //      mapping. To be investigated.
    const bool oldMoving = mesh_.moving(false);

    fvMeshAdder::add
    (
        Pstream::myProcNo(),    // index of mesh to modify (== mesh_)
        fvMeshes,
        oldFaceOwner,

        // Coupling info
        localBoundaryFace,
        remoteFaceProc,
        remoteBoundaryFace,

        constructPatchMap,
        constructCellMap,
        constructFaceMap,
        constructPointMap
    );

    mesh_.moving(oldMoving);


    if (debug)
    {
        Pout<< nl << "ADDED REMOTE MESHES:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printIntFieldInfo<volScalarField::Internal>(mesh_);
        printIntFieldInfo<volVectorField::Internal>(mesh_);
        printIntFieldInfo<volSphericalTensorField::Internal>(mesh_);
        printIntFieldInfo<volSymmTensorField::Internal>(mesh_);
        printIntFieldInfo<volTensorField::Internal>(mesh_);
        Pout<< nl << endl;
    }

    {
        //- Combine sourceProc, sourcePatch, sourceFace
        sourceProc.resize_nocopy(mesh_.nBoundaryFaces());
        sourceProc = -1;
        sourcePatch.resize_nocopy(mesh_.nBoundaryFaces());
        sourcePatch = -1;
        sourceFace.resize_nocopy(mesh_.nBoundaryFaces());
        sourceFace = -1;
        sourceNewNbrProc.resize_nocopy(mesh_.nBoundaryFaces());
        sourceNewNbrProc = -1;
        sourcePointMaster.resize_nocopy(mesh_.nPoints());
        sourcePointMaster = -1;

        if (mesh_.nPoints() > 0)
        {
            forAll(meshes, meshi)
            {
                if (domainSourceFaces.set(meshi))
                {
                    const label nIntFaces =
                    (
                        meshi == Pstream::myProcNo()
                      ? nOldInternalFaces
                      : meshes[meshi].nInternalFaces()
                    );
                    const labelList& faceOwner
                    (
                        meshi == Pstream::myProcNo()
                      ? oldFaceOwner
                      : meshes[meshi].faceOwner()
                    );

                    labelList& faceMap = constructFaceMap[meshi];
                    const labelList& cellMap = constructCellMap[meshi];

                    const labelList& domainSourceFace =
                        domainSourceFaces[meshi];
                    const labelList& domainSourceProc =
                        domainSourceProcs[meshi];
                    const labelList& domainSourcePatch =
                        domainSourcePatchs[meshi];
                    const labelList& domainSourceNewNbr =
                        domainSourceNewNbrProcs[meshi];
                    UIndirectList<label>
                    (
                        sourcePointMaster,
                        constructPointMap[meshi]
                    ) = domainSourcePointMasters[meshi];


                    forAll(domainSourceFace, bFacei)
                    {
                        const label oldFacei = bFacei+nIntFaces;
                        const label allFacei = faceMap[oldFacei];
                        const label allbFacei = allFacei-mesh_.nInternalFaces();

                        if (allbFacei >= 0)
                        {
                            sourceProc[allbFacei] = domainSourceProc[bFacei];
                            sourcePatch[allbFacei] = domainSourcePatch[bFacei];
                            sourceFace[allbFacei] = domainSourceFace[bFacei];
                            sourceNewNbrProc[allbFacei] =
                                domainSourceNewNbr[bFacei];
                        }
                    }


                    // Add flip to constructFaceMap
                    forAll(faceMap, oldFacei)
                    {
                        const label allFacei = faceMap[oldFacei];
                        const label allOwn = mesh_.faceOwner()[allFacei];

                        if (cellMap[faceOwner[oldFacei]] == allOwn)
                        {
                            // Master face
                            faceMap[oldFacei] += 1;
                        }
                        else
                        {
                            // Slave face. Flip.
                            faceMap[oldFacei] = -faceMap[oldFacei] - 1;
                        }
                    }
                }
            }
        }
    }


    UPstream::parRun(oldParRun);  // Restore parallel state


    // Print a bit.
    if (debug)
    {
        Pout<< nl << "REDISTRIBUTED MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printIntFieldInfo<volScalarField::Internal>(mesh_);
        printIntFieldInfo<volVectorField::Internal>(mesh_);
        printIntFieldInfo<volSphericalTensorField::Internal>(mesh_);
        printIntFieldInfo<volSymmTensorField::Internal>(mesh_);
        printIntFieldInfo<volTensorField::Internal>(mesh_);
        Pout<< nl << endl;
    }


    // See if any originally shared points need to be merged. Note: does
    // parallel comms. After this points and edges should again be consistent.
    mergeSharedPoints(sourcePointMaster, constructPointMap);


    // Add processorPatches
    // ~~~~~~~~~~~~~~~~~~~~

    // Per neighbour processor, per originating patch, the patchID
    // For faces resulting from internal faces or normal processor patches
    // the originating patch is -1. For cyclics this is the cyclic patchID.
    List<Map<label>> procPatchID;

    // Add processor and processorCyclic patches.
    addProcPatches(sourceNewNbrProc, sourcePatch, procPatchID);

    // Put faces into correct patch. Note that we now have proper
    // processorPolyPatches again so repatching will take care of coupled face
    // ordering.

    // Get boundary faces to be repatched. Is -1 or new patchID
    labelList newPatchID
    (
        getBoundaryPatch
        (
            sourceNewNbrProc,
            sourcePatch,
            procPatchID
        )
    );

    // Change patches. Since this might change ordering of coupled faces
    // we also need to adapt our constructMaps.
    repatch(newPatchID, constructFaceMap);

    // Bit of hack: processorFvPatchField does not get reset since created
    // from nothing so explicitly reset.
    initPatchFields<volScalarField, processorFvPatchField<scalar>>
    (
        Zero
    );
    initPatchFields<volVectorField, processorFvPatchField<vector>>
    (
        Zero
    );
    initPatchFields
    <
        volSphericalTensorField,
        processorFvPatchField<sphericalTensor>
    >
    (
        Zero
    );
    initPatchFields<volSymmTensorField, processorFvPatchField<symmTensor>>
    (
        Zero
    );
    initPatchFields<volTensorField, processorFvPatchField<tensor>>
    (
        Zero
    );


    mesh_.setInstance(mesh_.time().timeName());


    // Print a bit
    if (debug)
    {
        Pout<< nl << "FINAL MESH:" << endl;
        printMeshInfo(mesh_);
        printFieldInfo<volScalarField>(mesh_);
        printFieldInfo<volVectorField>(mesh_);
        printFieldInfo<volSphericalTensorField>(mesh_);
        printFieldInfo<volSymmTensorField>(mesh_);
        printFieldInfo<volTensorField>(mesh_);
        printFieldInfo<surfaceScalarField>(mesh_);
        printFieldInfo<surfaceVectorField>(mesh_);
        printFieldInfo<surfaceSphericalTensorField>(mesh_);
        printFieldInfo<surfaceSymmTensorField>(mesh_);
        printFieldInfo<surfaceTensorField>(mesh_);
        printIntFieldInfo<volScalarField::Internal>(mesh_);
        printIntFieldInfo<volVectorField::Internal>(mesh_);
        printIntFieldInfo<volSphericalTensorField::Internal>(mesh_);
        printIntFieldInfo<volSymmTensorField::Internal>(mesh_);
        printIntFieldInfo<volTensorField::Internal>(mesh_);
        Pout<< nl << endl;
    }

    // Collect all maps and return
    return autoPtr<mapDistributePolyMesh>::New
    (
        mesh_,

        nOldPoints,
        nOldFaces,
        nOldCells,
        std::move(oldPatchStarts),
        std::move(oldPatchNMeshPoints),

        std::move(subPointMap),
        std::move(subFaceMap),
        std::move(subCellMap),
        std::move(subPatchMap),

        std::move(constructPointMap),
        std::move(constructFaceMap),
        std::move(constructCellMap),
        std::move(constructPatchMap),

        true,           // subFaceMap has flip
        true            // constructFaceMap has flip
    );
}


// ************************************************************************* //
