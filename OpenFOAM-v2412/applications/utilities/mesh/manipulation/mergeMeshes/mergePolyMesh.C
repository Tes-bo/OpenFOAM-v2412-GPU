/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "mergePolyMesh.H"
#include "Time.H"
#include "polyTopoChanger.H"
#include "mapPolyMesh.H"
#include "polyAddPoint.H"
#include "polyAddCell.H"
#include "polyAddFace.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mergePolyMesh, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::mergePolyMesh::patchIndex(const polyPatch& p)
{
    // Find the patch name on the list.  If the patch is already there
    // and patch types match, return index
    const word& pType = p.type();
    const word& pName = p.name();

    bool nameFound = false;

    forAll(patchNames_, patchi)
    {
        if (patchNames_[patchi] == pName)
        {
            if (patchDicts_[patchi].get<word>("type") == pType)
            {
                // Found name and types match
                return patchi;
            }
            else
            {
                // Found the name, but type is different
                nameFound = true;
            }
        }
    }

    // Patch not found.  Append to the list
    {
        OCharStream os;
        p.write(os);
        ISpanStream is(os.view());

        patchDicts_.push_back(dictionary(is));
    }

    if (nameFound)
    {
        // Duplicate name is not allowed.  Create a composite name from the
        // patch name and case name
        const word& caseName = p.boundaryMesh().mesh().time().caseName();

        patchNames_.append(pName + "_" + caseName);

        Info<< "label patchIndex(const polyPatch& p) : "
            << "Patch " << p.index() << " named "
            << pName << " in mesh " << caseName
            << " already exists, but patch types "
            << " do not match.\nCreating a composite name as "
            << patchNames_.last() << endl;
    }
    else
    {
        patchNames_.append(pName);
    }

    return patchNames_.size() - 1;
}


Foam::label Foam::mergePolyMesh::zoneIndex
(
    DynamicList<word>& names,
    const word& curName
)
{
    forAll(names, zonei)
    {
        if (names[zonei] == curName)
        {
            return zonei;
        }
    }

    // Not found.  Add new name to the list
    names.append(curName);

    return names.size() - 1;
}


void Foam::mergePolyMesh::sortProcessorPatches()
{
    Info<< "Reordering processor patches last" << endl;

    // Updates boundaryMesh() and meshMod_ to guarantee processor patches
    // are last. This could be done inside the merge() but it is far easier
    // to do separately.


    // 1. Shuffle the patches in the boundaryMesh

    const polyBoundaryMesh& oldPatches = boundaryMesh();

    polyPatchList newPatches(oldPatches.size());

    labelList oldToSorted(oldPatches.size());

    label nPatches = 0;

    forAll(oldPatches, patchi)
    {
        const polyPatch& pp = oldPatches[patchi];

        if (!isA<processorPolyPatch>(pp))
        {
            newPatches.set
            (
                nPatches,
                pp.clone
                (
                    oldPatches,
                    nPatches,
                    0,
                    nInternalFaces()
                )
            );

            oldToSorted[patchi] = nPatches;
            ++nPatches;
        }
    }

    forAll(oldPatches, patchi)
    {
        const polyPatch& pp = oldPatches[patchi];

        if (isA<processorPolyPatch>(pp))
        {
            newPatches.set
            (
                nPatches,
                pp.clone
                (
                    oldPatches,
                    oldToSorted[patchi],
                    0,
                    nInternalFaces()
                )
            );

            oldToSorted[patchi] = nPatches;
            ++nPatches;
        }
    }


    removeBoundary();
    addPatches(newPatches);


    // Update the polyTopoChange
    DynamicList<label>& patchID = const_cast<DynamicList<label>&>
    (
        meshMod_.region()
    );

    forAll(patchID, facei)
    {
        label patchi = patchID[facei];
        if (patchi != -1)
        {
            patchID[facei] = oldToSorted[patchID[facei]];
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mergePolyMesh::mergePolyMesh(const IOobject& io)
:
    polyMesh(io),
    meshMod_(*this),
    patchNames_(2*boundaryMesh().size()),
    patchDicts_(2*boundaryMesh().size()),
    pointZoneNames_(),
    faceZoneNames_(),
    cellZoneNames_()
{
    // Insert the original patches into the list
    wordList curPatchNames = boundaryMesh().names();

    OCharStream os;
    forAll(boundaryMesh(), patchi)
    {
        os.rewind();
        boundaryMesh()[patchi].write(os);
        ISpanStream is(os.view());

        patchNames_.push_back(boundaryMesh()[patchi].name());
        patchDicts_.push_back(dictionary(is));
    }

    // Insert point, face and cell zones into the list

    // Point zones
    wordList curPointZoneNames = pointZones().names();
    if (curPointZoneNames.size())
    {
        pointZoneNames_.setCapacity(2*curPointZoneNames.size());
        pointZoneNames_.append(curPointZoneNames);
    }

    // Face zones
    wordList curFaceZoneNames = faceZones().names();
    if (curFaceZoneNames.size())
    {
        faceZoneNames_.setCapacity(2*curFaceZoneNames.size());
        faceZoneNames_.append(curFaceZoneNames);
    }

    // Cell zones
    wordList curCellZoneNames = cellZones().names();
    if (curCellZoneNames.size())
    {
        cellZoneNames_.setCapacity(2*curCellZoneNames.size());
        cellZoneNames_.append(curCellZoneNames);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mergePolyMesh::addMesh(const polyMesh& m)
{
    // Add all the points, faces and cells of the new mesh

    // Add points

    label zoneID = -1;

    const pointField& p = m.points();
    labelList renumberPoints(p.size());

    const pointZoneMesh& pz = m.pointZones();
    labelList pointZoneIndices(pz.size());

    forAll(pz, zoneI)
    {
        pointZoneIndices[zoneI] = zoneIndex(pointZoneNames_, pz[zoneI].name());
    }

    forAll(p, pointi)
    {
        // Grab zone ID.  If a point is not in a zone, it will return -1
        zoneID = pz.whichZone(pointi);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = pointZoneIndices[zoneID];
        }

        renumberPoints[pointi] =
            meshMod_.setAction
            (
                polyAddPoint
                (
                    p[pointi],            // Point to add
                    -1,                   // Master point (straight addition)
                    zoneID,               // Zone for point
                    pointi < m.nPoints()  // Is in cell?
                )
            );
    }

    // Add cells

    const cellList& c = m.cells();
    labelList renumberCells(c.size());

    const cellZoneMesh& cz = m.cellZones();
    labelList cellZoneIndices(cz.size());

    forAll(cz, zoneI)
    {
        cellZoneIndices[zoneI] = zoneIndex(cellZoneNames_, cz[zoneI].name());
    }

    forAll(c, celli)
    {
        // Grab zone ID.  If a cell is not in a zone, it will return -1
        zoneID = cz.whichZone(celli);

        if (zoneID >= 0)
        {
            // Translate zone ID into the new index
            zoneID = cellZoneIndices[zoneID];
        }

        renumberCells[celli] =
            meshMod_.setAction
            (
                polyAddCell
                (
                    -1,                   // Master point
                    -1,                   // Master edge
                    -1,                   // Master face
                    -1,                   // Master cell
                    zoneID                // Zone for cell
                )
            );
    }

    // Add faces
    const polyBoundaryMesh& bm = m.boundaryMesh();

    // Gather the patch indices
    labelList patchIndices(bm.size());

    forAll(patchIndices, patchi)
    {
        patchIndices[patchi] = patchIndex(bm[patchi]);
    }

    // Temporary: update number of allowable patches. This should be
    // determined at the top - before adding anything.
    meshMod_.setNumPatches(patchNames_.size());



    const faceZoneMesh& fz = m.faceZones();
    labelList faceZoneIndices(fz.size());

    forAll(fz, zoneI)
    {
        faceZoneIndices[zoneI] = zoneIndex(faceZoneNames_, fz[zoneI].name());
    }

    const faceList& f = m.faces();
    labelList renumberFaces(f.size());

    const labelList& own = m.faceOwner();
    const labelList& nei = m.faceNeighbour();

    label newOwn, newNei, newPatch, newZone;
    bool newZoneFlip;

    forAll(f, facei)
    {
        const face& curFace = f[facei];

        face newFace(curFace.size());

        forAll(curFace, pointi)
        {
            newFace[pointi] = renumberPoints[curFace[pointi]];
        }

        if (debug)
        {
            // Check that the face is valid
            if (min(newFace) < 0)
            {
                FatalErrorInFunction
                    << "Error in point mapping for face " << facei
                    << ".  Old face: " << curFace << " New face: " << newFace
                    << abort(FatalError);
            }
        }

        if (facei < m.nInternalFaces() || facei >= m.nFaces())
        {
            newPatch = -1;
        }
        else
        {
            newPatch = patchIndices[bm.whichPatch(facei)];
        }

        newOwn = own[facei];
        if (newOwn > -1) newOwn = renumberCells[newOwn];

        if (newPatch > -1)
        {
            newNei = -1;
        }
        else
        {
            newNei = nei[facei];
            newNei = renumberCells[newNei];
        }


        newZone = fz.whichZone(facei);
        newZoneFlip = false;

        if (newZone >= 0)
        {
            newZoneFlip = fz[newZone].flipMap()[fz[newZone].whichFace(facei)];

            // Grab the new zone
            newZone = faceZoneIndices[newZone];
        }

        renumberFaces[facei] =
            meshMod_.setAction
            (
                polyAddFace
                (
                    newFace,
                    newOwn,
                    newNei,
                    -1,
                    -1,
                    -1,
                    false,
                    newPatch,
                    newZone,
                    newZoneFlip
                )
            );
    }
}


void Foam::mergePolyMesh::merge()
{
    Info<< "patch names: " << patchNames_ << nl
        << "patch dicts: " << patchDicts_ << nl
        << "point zone names: " << pointZoneNames_ << nl
        << "face zone names: " << faceZoneNames_ << nl
        << "cell zone names: " << cellZoneNames_ << endl;

    // Add the patches if necessary
    if (patchNames_.size() != boundaryMesh().size())
    {
        Info<< "Copying old patches" << endl;

        polyPatchList newPatches(patchNames_.size());

        const polyBoundaryMesh& oldPatches = boundaryMesh();

        // Note.  Re-using counter in two for loops
        label patchi = 0;

        for (patchi = 0; patchi < oldPatches.size(); patchi++)
        {
            newPatches.set
            (
                patchi,
                oldPatches[patchi].clone(oldPatches)
            );
        }

        Info<< "Adding new patches. " << endl;

        label endOfLastPatch =
            patchi == 0
          ? 0
          : oldPatches[patchi - 1].start() + oldPatches[patchi - 1].size();

        for (; patchi < patchNames_.size(); patchi++)
        {
            // Add a patch
            dictionary dict(patchDicts_[patchi]);
            dict.set("nFaces", 0);
            dict.set("startFace", endOfLastPatch);

            newPatches.set
            (
                patchi,
                (
                    polyPatch::New
                    (
                        patchNames_[patchi],
                        dict,
                        patchi,
                        oldPatches
                    )
                )
            );
        }

        removeBoundary();
        addPatches(newPatches);
    }

    // Add the zones if necessary
    if (pointZoneNames_.size() > pointZones().size())
    {
        Info<< "Adding new pointZones." << endl;

        label zonei = pointZones().size(); // continue from here

        const label nZones = pointZoneNames_.size();

        pointZones().setSize(nZones);

        for (/*nil*/; zonei < nZones; ++zonei)
        {
            pointZones().set
            (
                zonei,
                new pointZone(pointZoneNames_[zonei], zonei, pointZones())
            );
        }
    }
    if (cellZoneNames_.size() > cellZones().size())
    {
        Info<< "Adding new cellZones." << endl;

        label zonei = cellZones().size(); // continue from here

        const label nZones = cellZoneNames_.size();

        cellZones().setSize(cellZoneNames_.size());

        for (/*nil*/; zonei < nZones; ++zonei)
        {
            cellZones().set
            (
                zonei,
                new cellZone(cellZoneNames_[zonei], zonei, cellZones())
            );
        }
    }
    if (faceZoneNames_.size() > faceZones().size())
    {
        Info<< "Adding new faceZones." << endl;

        label zonei = faceZones().size(); // continue from here

        const label nZones = faceZoneNames_.size();

        faceZones().setSize(nZones);

        for (/*nil*/; zonei < nZones; ++zonei)
        {
            faceZones().set
            (
                zonei,
                new faceZone(faceZoneNames_[zonei], zonei, faceZones())
            );
        }
    }


    // Shuffle the processor patches to be last
    sortProcessorPatches();

    // Change mesh. No inflation
    meshMod_.changeMesh(*this, false);

    // Clear topo change for the next operation
    meshMod_.clear();
}


// ************************************************************************* //
