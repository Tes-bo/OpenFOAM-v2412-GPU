/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "slidingInterface.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"
#include "polyTopoChanger.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::calcAttachedAddressing() const
{
    if (debug)
    {
        Pout<< FUNCTION_NAME
            << " for object " << name() << " : "
            << "Calculating zone face-cell addressing."
            << endl;
    }

    if (!attached_)
    {
        // Clear existing addressing
        clearAttachedAddressing();

        const polyMesh& mesh = topoChanger().mesh();
        const labelList& own = mesh.faceOwner();
        const labelList& nei = mesh.faceNeighbour();
        const faceZoneMesh& faceZones = mesh.faceZones();

        // Master side

        const primitiveFacePatch& masterPatch =
            faceZones[masterFaceZoneID_.index()]();

        const labelList& masterPatchFaces =
            faceZones[masterFaceZoneID_.index()];

        const boolList& masterFlip =
            faceZones[masterFaceZoneID_.index()].flipMap();

        masterFaceCellsPtr_.reset(new labelList(masterPatchFaces.size()));
        auto& mfc = *masterFaceCellsPtr_;

        forAll(masterPatchFaces, facei)
        {
            if (masterFlip[facei])
            {
                mfc[facei] = nei[masterPatchFaces[facei]];
            }
            else
            {
                mfc[facei] = own[masterPatchFaces[facei]];
            }
        }

        // Slave side

        const primitiveFacePatch& slavePatch =
            faceZones[slaveFaceZoneID_.index()]();

        const labelList& slavePatchFaces =
            faceZones[slaveFaceZoneID_.index()];

        const boolList& slaveFlip =
            faceZones[slaveFaceZoneID_.index()].flipMap();

        slaveFaceCellsPtr_.reset(new labelList(slavePatchFaces.size()));
        auto& sfc = *slaveFaceCellsPtr_;

        forAll(slavePatchFaces, facei)
        {
            if (slaveFlip[facei])
            {
                sfc[facei] = nei[slavePatchFaces[facei]];
            }
            else
            {
                sfc[facei] = own[slavePatchFaces[facei]];
            }
        }

        // Check that the addressing is valid
        if (min(mfc) < 0 || min(sfc) < 0)
        {
            if (debug)
            {
                forAll(mfc, facei)
                {
                    if (mfc[facei] < 0)
                    {
                        Pout<< "No cell next to master patch face " << facei
                            << ".  Global face no: " << mfc[facei]
                            << " own: " << own[masterPatchFaces[facei]]
                            << " nei: " << nei[masterPatchFaces[facei]]
                            << " flip: " << masterFlip[facei] << endl;
                    }
                }

                forAll(sfc, facei)
                {
                    if (sfc[facei] < 0)
                    {
                        Pout<< "No cell next to slave patch face " << facei
                            << ".  Global face no: " << sfc[facei]
                            << " own: " << own[slavePatchFaces[facei]]
                            << " nei: " << nei[slavePatchFaces[facei]]
                            << " flip: " << slaveFlip[facei] << endl;
                    }
                }
            }

            FatalErrorInFunction
                << "decoupled mesh or sliding interface definition."
                << abort(FatalError);
        }

        // Calculate stick-out faces
        const labelListList& pointFaces = mesh.pointFaces();

        labelHashSet stickOutFaceMap
        (
            primitiveMesh::facesPerCell_
          * max(masterPatch.size(), slavePatch.size())
        );

        // Master side
        const labelList& masterMeshPoints = masterPatch.meshPoints();

        stickOutFaceMap.clear();

        for (const label pointi : masterMeshPoints)
        {
            for (const label facei : pointFaces[pointi])
            {
                const label zoneIdx = faceZones.whichZone(facei);

                // Add if face not already part of master or slave face zone
                // This handles partially attached faces.
                if
                (
                    zoneIdx != masterFaceZoneID_.index()
                 && zoneIdx != slaveFaceZoneID_.index()
                )
                {
                    stickOutFaceMap.insert(facei);
                }
            }
        }

        masterStickOutFacesPtr_.reset(new labelList(stickOutFaceMap.toc()));

        // Sort in debug mode for easier diagnostics
        if (debug)
        {
            Foam::sort(*masterStickOutFacesPtr_);
        }

        // Slave side
        const labelList& slaveMeshPoints = slavePatch.meshPoints();

        stickOutFaceMap.clear();

        for (const label pointi : slaveMeshPoints)
        {
            for (const label facei : pointFaces[pointi])
            {
                const label zoneIdx = faceZones.whichZone(facei);

                // Add if face not already part of master or slave face zone
                // This handles partially attached faces.
                if
                (
                    zoneIdx != masterFaceZoneID_.index()
                 && zoneIdx != slaveFaceZoneID_.index()
                )
                {
                    stickOutFaceMap.insert(facei);
                }
            }
        }

        slaveStickOutFacesPtr_.reset(new labelList(stickOutFaceMap.toc()));

        // Sort in debug mode for easier diagnostics
        if (debug)
        {
            Foam::sort(*slaveStickOutFacesPtr_);
        }

        stickOutFaceMap.clear();

        // Retired point addressing does not exist at this stage.
        // It will be filled when the interface is coupled.
        retiredPointMapPtr_.reset
        (
            new Map<label>
            (
                2*faceZones[slaveFaceZoneID_.index()]().nPoints()
            )
        );

        // Ditto for cut point edge map.  This is a rough guess of its size
        cutPointEdgePairMapPtr_.reset
        (
            new Map<Pair<edge>>
            (
                faceZones[slaveFaceZoneID_.index()]().nEdges()
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "cannot be assembled for object " << name()
            << abort(FatalError);
    }

    if (debug)
    {
        Pout<< FUNCTION_NAME
            << " for object " << name() << " : "
            << "Finished calculating zone face-cell addressing."
            << endl;
    }
}


void Foam::slidingInterface::clearAttachedAddressing() const
{
    masterFaceCellsPtr_.reset(nullptr);
    slaveFaceCellsPtr_.reset(nullptr);

    masterStickOutFacesPtr_.reset(nullptr);
    slaveStickOutFacesPtr_.reset(nullptr);

    retiredPointMapPtr_.reset(nullptr);
    cutPointEdgePairMapPtr_.reset(nullptr);
}


void Foam::slidingInterface::renumberAttachedAddressing
(
    const mapPolyMesh& m
) const
{
    // Get reference to reverse cell renumbering
    // The renumbering map is needed the other way around, i.e. giving
    // the new cell number for every old cell next to the interface.
    const labelList& reverseCellMap = m.reverseCellMap();

    const labelList& mfc = masterFaceCells();
    const labelList& sfc = slaveFaceCells();

    // Master side
    std::unique_ptr<labelList> newMfcPtr(new labelList(mfc.size(), -1));
    auto& newMfc = *newMfcPtr;

    const labelList& mfzRenumber =
        m.faceZoneFaceMap()[masterFaceZoneID_.index()];

    forAll(mfc, facei)
    {
        label newCelli = reverseCellMap[mfc[mfzRenumber[facei]]];

        if (newCelli >= 0)
        {
            newMfc[facei] = newCelli;
        }
    }

    // Slave side
    std::unique_ptr<labelList> newSfcPtr(new labelList(sfc.size(), -1));
    auto& newSfc = *newSfcPtr;

    const labelList& sfzRenumber =
        m.faceZoneFaceMap()[slaveFaceZoneID_.index()];

    forAll(sfc, facei)
    {
        label newCelli = reverseCellMap[sfc[sfzRenumber[facei]]];

        if (newCelli >= 0)
        {
            newSfc[facei] = newCelli;
        }
    }

    if (debug)
    {
        // Check if all the mapped cells are live
        if (min(newMfc) < 0 || min(newSfc) < 0)
        {
            FatalErrorInFunction
                << "Error in cell renumbering for object " << name()
                << ".  Some of master cells next "
                << "to the interface have been removed."
                << abort(FatalError);
        }
    }

    // Renumber stick-out faces

    const labelList& reverseFaceMap = m.reverseFaceMap();

    // Master side
    const labelList& msof = masterStickOutFaces();

    std::unique_ptr<labelList> newMsofPtr(new labelList(msof.size(), -1));
    auto& newMsof = *newMsofPtr;

    forAll(msof, facei)
    {
        label newFacei = reverseFaceMap[msof[facei]];

        if (newFacei >= 0)
        {
            newMsof[facei] = newFacei;
        }
    }
//     Pout<< "newMsof: " << newMsof << endl;
    // Slave side
    const labelList& ssof = slaveStickOutFaces();

    std::unique_ptr<labelList> newSsofPtr(new labelList(ssof.size(), -1));
    auto& newSsof = *newSsofPtr;

    forAll(ssof, facei)
    {
        label newFacei = reverseFaceMap[ssof[facei]];

        if (newFacei >= 0)
        {
            newSsof[facei] = newFacei;
        }
    }
//     Pout<< "newSsof: " << newSsof << endl;
    if (debug)
    {
        // Check if all the mapped cells are live
        if (min(newMsof) < 0 || min(newSsof) < 0)
        {
            FatalErrorInFunction
                << "Error in face renumbering for object " << name()
                << ".  Some of stick-out next "
                << "to the interface have been removed."
                << abort(FatalError);
        }
    }

    // Renumber the retired point map. Need to take a copy!
    const Map<label> rpm = retiredPointMap();

    std::unique_ptr<Map<label>> newRpmPtr(new Map<label>(rpm.size()));
    auto& newRpm = *newRpmPtr;

    // Get reference to point renumbering
    const labelList& reversePointMap = m.reversePointMap();

    forAllConstIters(rpm, iter)
    {
        const label key = reversePointMap[iter.key()];
        const label val = reversePointMap[iter.val()];

        if (debug)
        {
            // Check if all the mapped cells are live
            if (key < 0 || val < 0)
            {
                FatalErrorInFunction
                    << "Error in retired point numbering for object "
                    << name() << ".  Some of master "
                    << "points have been removed."
                    << abort(FatalError);
            }
        }

        newRpm.insert(key, val);
    }

    // Renumber the cut point edge pair map. Need to take a copy!
    const Map<Pair<edge>> cpepm = cutPointEdgePairMap();

    std::unique_ptr<Map<Pair<edge>>> newCpepmPtr
    (
        new Map<Pair<edge>>(cpepm.size())
    );
    auto& newCpepm = *newCpepmPtr;

    forAllConstIters(cpepm, iter)
    {
        const label key = reversePointMap[iter.key()];

        const Pair<edge>& oldPe = iter.val();

        // Re-do the edges in global addressing
        const label ms = reversePointMap[oldPe.first().start()];
        const label me = reversePointMap[oldPe.first().end()];

        const label ss = reversePointMap[oldPe.second().start()];
        const label se = reversePointMap[oldPe.second().end()];

        if (debug)
        {
            // Check if all the mapped cells are live
            if (key < 0 || ms < 0 || me < 0 || ss < 0 || se < 0)
            {
                FatalErrorInFunction
                    << "Error in cut point edge pair map numbering for object "
                    << name() << ".  Some of master points have been removed."
                    << abort(FatalError);
            }
        }

        newCpepm.insert(key, Pair<edge>(edge(ms, me), edge(ss, se)));
    }

    if (!projectedSlavePointsPtr_)
    {
        FatalErrorInFunction
            << "Error in projected point numbering for object " << name()
            << abort(FatalError);
    }

    // Renumber the projected slave zone points
    const pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

    std::unique_ptr<pointField> newProjectedSlavePointsPtr
    (
        new pointField(projectedSlavePoints.size())
    );
    auto& newProjectedSlavePoints = *newProjectedSlavePointsPtr;

    const labelList& sfzPointRenumber =
        m.faceZonePointMap()[slaveFaceZoneID_.index()];

    forAll(newProjectedSlavePoints, pointi)
    {
        if (sfzPointRenumber[pointi] > -1)
        {
            newProjectedSlavePoints[pointi] =
                projectedSlavePoints[sfzPointRenumber[pointi]];
        }
    }

    // Re-set the lists
    clearAttachedAddressing();

    projectedSlavePointsPtr_.reset(nullptr);

    masterFaceCellsPtr_ = std::move(newMfcPtr);
    slaveFaceCellsPtr_ = std::move(newSfcPtr);

    masterStickOutFacesPtr_ = std::move(newMsofPtr);
    slaveStickOutFacesPtr_ = std::move(newSsofPtr);

    retiredPointMapPtr_ = std::move(newRpmPtr);
    cutPointEdgePairMapPtr_ = std::move(newCpepmPtr);
    projectedSlavePointsPtr_ = std::move(newProjectedSlavePointsPtr);
}


const Foam::labelList& Foam::slidingInterface::masterFaceCells() const
{
    if (!masterFaceCellsPtr_)
    {
        FatalErrorInFunction
            << "Master zone face-cell addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *masterFaceCellsPtr_;
}


const Foam::labelList& Foam::slidingInterface::slaveFaceCells() const
{
    if (!slaveFaceCellsPtr_)
    {
        FatalErrorInFunction
            << "Slave zone face-cell addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *slaveFaceCellsPtr_;
}


const Foam::labelList& Foam::slidingInterface::masterStickOutFaces() const
{
    if (!masterStickOutFacesPtr_)
    {
        FatalErrorInFunction
            << "Master zone stick-out face addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *masterStickOutFacesPtr_;
}


const Foam::labelList& Foam::slidingInterface::slaveStickOutFaces() const
{
    if (!slaveStickOutFacesPtr_)
    {
        FatalErrorInFunction
            << "Slave zone stick-out face addressing not available for object "
            << name()
            << abort(FatalError);
    }

    return *slaveStickOutFacesPtr_;
}


const Foam::Map<Foam::label>& Foam::slidingInterface::retiredPointMap() const
{
    if (!retiredPointMapPtr_)
    {
        FatalErrorInFunction
            << "Retired point map not available for object " << name()
            << abort(FatalError);
    }

    return *retiredPointMapPtr_;
}


const Foam::Map<Foam::Pair<Foam::edge>>&
Foam::slidingInterface::cutPointEdgePairMap() const
{
    if (!cutPointEdgePairMapPtr_)
    {
        FatalErrorInFunction
            << "Retired point map not available for object " << name()
            << abort(FatalError);
    }

    return *cutPointEdgePairMapPtr_;
}


// ************************************************************************* //
