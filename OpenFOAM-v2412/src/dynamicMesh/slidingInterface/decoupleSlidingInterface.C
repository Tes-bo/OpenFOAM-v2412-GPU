/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "polyTopoChanger.H"
#include "polyModifyFace.H"
#include "polyModifyPoint.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::decoupleInterface
(
    polyTopoChange& ref
) const
{
    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << ": Decoupling sliding interface " << name() << endl;
    }

    if (!attached_)
    {
        if (debug)
        {
            Pout<< FUNCTION_NAME << nl
                << ": Interface already decoupled." << endl;
        }

        return;
    }

    // Clear previous couple
    clearCouple(ref);

    const polyMesh& mesh = topoChanger().mesh();
    const pointField& points = mesh.points();
    const faceList& faces = mesh.faces();
    const cellList& cells = mesh.cells();
    const labelList& own = mesh.faceOwner();
    const labelList& nei = mesh.faceNeighbour();
    const faceZoneMesh& faceZones = mesh.faceZones();

    // Master side

    const primitiveFacePatch& masterPatch =
        faceZones[masterFaceZoneID_.index()]();

    const labelList& masterPatchAddr =
        faceZones[masterFaceZoneID_.index()];

    const boolList& masterPatchFlip =
        faceZones[masterFaceZoneID_.index()].flipMap();

    const labelList& masterFc = masterFaceCells();

    // Recover faces in master patch

    forAll(masterPatchAddr, facei)
    {
        // Make a copy of the face and turn it if necessary
        face newFace = faces[masterPatchAddr[facei]];

        if (masterPatchFlip[facei])
        {
            newFace.flip();
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,                         // new face
                masterPatchAddr[facei],          // master face index
                masterFc[facei],                 // owner
                -1,                              // neighbour
                false,                           // flux flip
                masterPatchID_.index(),          // patch ID
                false,                           // remove from zone
                masterFaceZoneID_.index(),       // zone ID
                false                            // zone flip.  Face corrected
            )
        );

        // Pout<< "Modifying master patch face no "
        //     << masterPatchAddr[facei]
        //     << " face: " << faces[masterPatchAddr[facei]]
        //     << " old owner: " << own[masterPatchAddr[facei]]
        //     << " new owner: " << masterFc[facei]
        //     << endl;
    }

    // Slave side

    const primitiveFacePatch& slavePatch =
        faceZones[slaveFaceZoneID_.index()]();

    const labelList& slavePatchAddr =
        faceZones[slaveFaceZoneID_.index()];

    const boolList& slavePatchFlip =
        faceZones[slaveFaceZoneID_.index()].flipMap();

    const labelList& slaveFc = slaveFaceCells();

    // Grab retired point mapping
    const Map<label>& rpm = retiredPointMap();

    // Recover faces in slave patch

    forAll(slavePatchAddr, facei)
    {
        // Make a copy of face and turn it if necessary
        face newFace = faces[slavePatchAddr[facei]];

        if (slavePatchFlip[facei])
        {
            newFace.flip();
        }

        // Recover retired points on the slave side
        forAll(newFace, pointi)
        {
            newFace[pointi] = rpm.lookup(newFace[pointi], newFace[pointi]);
        }

        ref.setAction
        (
            polyModifyFace
            (
                newFace,                         // new face
                slavePatchAddr[facei],           // master face index
                slaveFc[facei],                  // owner
                -1,                              // neighbour
                false,                           // flux flip
                slavePatchID_.index(),           // patch ID
                false,                           // remove from zone
                slaveFaceZoneID_.index(),        // zone ID
                false                            // zone flip.  Face corrected
            )
        );
    }

    // Re-create the master stick-out faces

    // Grab the list of faces in the layer
    const labelList& masterStickOuts = masterStickOutFaces();

    for (const label curFaceID : masterStickOuts)
    {
        // Renumber the face and remove additional points

        const face& oldFace = faces[curFaceID];

        DynamicList<label> newFaceLabels(oldFace.size());

        bool changed = false;

        forAll(oldFace, pointi)
        {
            // Check if the point is removed
            if (ref.pointRemoved(oldFace[pointi]))
            {
                // Point removed; skip it
                changed = true;
            }
            else
            {
                newFaceLabels.append(oldFace[pointi]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorInFunction
                    << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            const label modifiedFaceZone = faceZones.whichZone(curFaceID);

            const bool modifiedFaceZoneFlip =
            (
                modifiedFaceZone >= 0
              ?
                faceZones[modifiedFaceZone].flipMap()
                [
                    faceZones[modifiedFaceZone].whichFace(curFaceID)
                ]
              : false
            );

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying master stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    curFaceID,              // label of face being modified
                    own[curFaceID],         // owner
                    nei[curFaceID],         // neighbour
                    false,                  // face flip
                    mesh.boundaryMesh().whichPatch(curFaceID), // patch for face
                    false,                  // remove from zone
                    modifiedFaceZone,       // zone for face
                    modifiedFaceZoneFlip    // face flip in zone
                )
            );
        }
    }

    // Re-create the slave stick-out faces

    labelHashSet slaveLayerCellFaceMap
    (
        primitiveMesh::facesPerCell_*(masterPatch.size() + slavePatch.size())
    );

    for (const label slaveFci : slaveFc)
    {
        const labelList& curFaces = cells[slaveFci];

        for (const label facei : curFaces)
        {
            // Check if the face belongs to the slave face zone; and
            // if it has been removed; if not add it
            if
            (
                faceZones.whichZone(facei)
             != slaveFaceZoneID_.index()
             && !ref.faceRemoved(facei)
            )
            {
                slaveLayerCellFaceMap.insert(facei);
            }
        }
    }

    // Grab the list of faces in the layer
    const labelList& slaveStickOuts = slaveStickOutFaces();

    // Grab master point mapping
    const Map<label>& masterPm = masterPatch.meshPointMap();

    for (const label curFaceID : slaveStickOuts)
    {
        // Renumber the face and remove additional points

        const face& oldFace = faces[curFaceID];

        DynamicList<label> newFaceLabels(oldFace.size());

        bool changed = false;

        forAll(oldFace, pointi)
        {
            // Check if the point is removed or retired

            const label retiredPointi = rpm.lookup(oldFace[pointi], -1);

            if (retiredPointi != -1)
            {
                // Master of retired point; grab its original
                changed = true;

                // Pout<< "Reinstating retired point: " << oldFace[pointi]
                //     << " with old: " << retiredPointi
                //     << endl;

                newFaceLabels.append(retiredPointi);
            }
            else if (ref.pointRemoved(oldFace[pointi]))
            {
                // Point removed; skip it
                changed = true;
            }
            else if (masterPm.found(oldFace[pointi]))
            {
                // Point from master patch only; skip it
                changed = true;
            }
            else
            {
                newFaceLabels.append(oldFace[pointi]);
            }
        }

        if (changed)
        {
            if (newFaceLabels.size() < 3)
            {
                FatalErrorInFunction
                    << "Face " << curFaceID << " reduced to less than "
                    << "3 points.  Topological/cutting error." << nl
                    << "Old face: " << oldFace << " new face: " << newFaceLabels
                    << abort(FatalError);
            }

            // Get face zone and its flip
            const label modifiedFaceZone =
                faceZones.whichZone(curFaceID);

            const bool modifiedFaceZoneFlip =
            (
                modifiedFaceZone >= 0
              ?
                faceZones[modifiedFaceZone].flipMap()
                [
                    faceZones[modifiedFaceZone].whichFace(curFaceID)
                ]
              : false
            );

            face newFace;
            newFace.transfer(newFaceLabels);

            // Pout<< "Modifying slave stick-out face " << curFaceID
            //     << " old face: " << oldFace
            //     << " new face: " << newFace
            //     << endl;

            // Modify the face
            ref.setAction
            (
                polyModifyFace
                (
                    newFace,                // modified face
                    curFaceID,              // label of face being modified
                    own[curFaceID],         // owner
                    nei[curFaceID],         // neighbour
                    false,                  // face flip
                    mesh.boundaryMesh().whichPatch(curFaceID), // patch for face
                    false,                  // remove from zone
                    modifiedFaceZone,       // zone for face
                    modifiedFaceZoneFlip    // face flip in zone
                )
            );
        }
    }

    // Bring all slave patch points back to life
    const labelList& slaveMeshPoints =
        faceZones[slaveFaceZoneID_.index()]().meshPoints();

    for (const label slavePointi : slaveMeshPoints)
    {
        ref.setAction
        (
            polyModifyPoint
            (
                slavePointi,            // point ID
                points[slavePointi],    // point
                false,                  // remove from zone
                mesh.pointZones().whichZone(slavePointi), // zone
                true                    // in a cell
            )
        );
    }

    // Clear the retired point numbering
    retiredPointMapPtr_->clear();

    // Finished decoupling
    attached_ = false;

    if (debug)
    {
        Pout<< FUNCTION_NAME << nl
            << ": Finished decoupling sliding interface " << name() << endl;
    }
}


// ************************************************************************* //
