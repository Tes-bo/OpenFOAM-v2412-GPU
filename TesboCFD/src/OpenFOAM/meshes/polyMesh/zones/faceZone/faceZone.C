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

#include "faceZone.H"
#include "addToRunTimeSelectionTable.H"
#include "faceZoneMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "mapPolyMesh.H"
#include "syncTools.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceZone, 0);
    defineRunTimeSelectionTable(faceZone, dictionary);
    addToRunTimeSelectionTable(faceZone, faceZone, dictionary);
}

const char* const Foam::faceZone::labelsName = "faceLabels";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceZone::setFlipMap(const bool val)
{
    // Match size for flipMap
    flipMap_.resize_nocopy(this->size());
    flipMap_ = val;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::faceZone::calcFaceZonePatch() const
{
    DebugInFunction << "Calculating primitive patch" << endl;

    if (patchPtr_)
    {
        FatalErrorInFunction
            << "primitive face zone patch already calculated"
            << abort(FatalError);
    }

    patchPtr_.reset
    (
        new primitiveFacePatch
        (
            faceList(size()),
            zoneMesh().mesh().points()
        )
    );
    auto& patch = *patchPtr_;

    const faceList& f = zoneMesh().mesh().faces();

    const labelList& addr = *this;
    const boolList& flips = flipMap();

    forAll(addr, facei)
    {
        if (flips[facei])
        {
            patch[facei] = f[addr[facei]].reverseFace();
        }
        else
        {
            patch[facei] = f[addr[facei]];
        }
    }

    DebugInfo << "Finished calculating primitive patch" << endl;
}


void Foam::faceZone::calcCellLayers() const
{
    DebugInFunction << "Calculating cell layers" << endl;

    if (frontCellsPtr_ || backCellsPtr_)
    {
        FatalErrorInFunction
            << "cell layers already calculated"
            << abort(FatalError);
    }
    else
    {
        // Go through all the faces in the zone.
        // Choose the front/back cell based on the face flip

        const labelList& own = zoneMesh().mesh().faceOwner();
        const labelList& nei = zoneMesh().mesh().faceNeighbour();

        const labelList& addr = *this;
        const boolList& flips = flipMap();

        frontCellsPtr_.reset(new labelList(addr.size()));
        backCellsPtr_.reset(new labelList(addr.size()));

        auto& fronts = *frontCellsPtr_;
        auto& backs = *backCellsPtr_;

        forAll(addr, facei)
        {
            const label ownCelli = own[addr[facei]];
            const label neiCelli =
            (
                zoneMesh().mesh().isInternalFace(addr[facei])
              ? nei[addr[facei]]
              : -1
            );

            if (flips[facei])
            {
                fronts[facei] = ownCelli;
                backs[facei] = neiCelli;
            }
            else
            {
                fronts[facei] = neiCelli;
                backs[facei] = ownCelli;
            }
        }
    }
}


// Foam::label Foam::faceZone::getLayerCell
// (
//     const side which,
//     const label i
// ) const
// {
//     const label facei = labelList::operator[](i);
//     const bool flipped = flipMap_[i];
//
//     if (which == side::FRONT ? flipped : !flipped)
//     {
//         return zoneMesh().mesh().faceOwner()[facei];
//     }
//     else if (zoneMesh().mesh().isInternalFace(facei))
//     {
//         return zoneMesh().mesh().faceNeighbour()[facei];
//     }
//     else
//     {
//         return -1;
//     }
// }
//
//
// Foam::label Foam::faceZone::frontCell(const label i) const
// {
//     return getLayerCell(side::FRONT, i);
// }
//
//
// Foam::label Foam::faceZone::backCell(const label i) const
// {
//     return getLayerCell(side::BACK, i);
// }


void Foam::faceZone::checkAddressing() const
{
    const labelList& addr = *this;

    if (addr.size() != flipMap_.size())
    {
        FatalErrorInFunction
            << "Size of addressing: " << addr.size()
            << " size of flip map: " << flipMap_.size()
            << abort(FatalError);
    }

    // Note: nFaces, nCells might not be set yet on mesh so use owner size
    const label nFaces = zoneMesh().mesh().faceOwner().size();

    for (const label facei : addr)
    {
        if (facei < 0 || facei >= nFaces)
        {
            WarningInFunction
                << "Illegal face index " << facei
                << " outside range 0.." << nFaces-1 << endl;
            break;  // Only report once
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faceZone::faceZone(const faceZoneMesh& zm)
:
    faceZone(word::null, 0, zm)
{}


Foam::faceZone::faceZone
(
    const word& name,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(name, index),
    zoneMesh_(zm)
{}


Foam::faceZone::faceZone
(
    const word& name,
    const labelUList& addr,
    const bool flipMapValue,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(name, index, zm)
{
    labelList::operator=(addr);
    flipMap_.resize(labelList::size(), flipMapValue);

    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    labelList&& addr,
    const bool flipMapValue,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(name, index, zm)
{
    labelList::transfer(addr);
    flipMap_.resize(labelList::size(), flipMapValue);
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const labelUList& addr,
    const boolUList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(name, index, zm)
{
    labelList::operator=(addr);
    flipMap_ = fm;

    // TBD
    // if (flipMap_.empty())
    // {
    //     // An empty flipMap is treated like 'false' instead of as an error
    //     flipMap_.resize(labelList::size(), false);
    // }

    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    labelList&& addr,
    boolList&& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(name, index, zm)
{
    labelList::transfer(addr);
    flipMap_.transfer(fm);

    // TBD
    // if (flipMap_.empty())
    // {
    //     // An empty flipMap is treated like 'false' instead of as an error
    //     flipMap_.resize(labelList::size(), false);
    // }

    checkAddressing();
}


Foam::faceZone::faceZone
(
    const word& name,
    const dictionary& dict,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(name, dict, this->labelsName, index),
    flipMap_(dict.lookup("flipMap")),  // OR: dict.get<boolList>("flipMap")
    zoneMesh_(zm)
{
    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& originalZone,
    const Foam::zero,
    const faceZoneMesh& zm,
    const label newIndex
)
:
    zone(originalZone, labelList(), newIndex),
    zoneMesh_(zm)
{}


Foam::faceZone::faceZone
(
    const faceZone& originalZone,
    const Foam::zero,
    const label index,
    const faceZoneMesh& zm
)
:
    zone(originalZone, labelList(), index),
    zoneMesh_(zm)
{}


Foam::faceZone::faceZone
(
    const faceZone& originalZone,
    const labelUList& addr,
    const boolUList& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(originalZone, Foam::zero{}, index, zm)
{
    labelList::operator=(addr);
    flipMap_ = fm;

    checkAddressing();
}


Foam::faceZone::faceZone
(
    const faceZone& originalZone,
    labelList&& addr,
    boolList&& fm,
    const label index,
    const faceZoneMesh& zm
)
:
    faceZone(originalZone, Foam::zero{}, index, zm)
{
    labelList::transfer(addr);
    flipMap_.transfer(fm);

    checkAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::faceZone::whichFace(const label globalFaceID) const
{
    return zone::localID(globalFaceID);
}


const Foam::primitiveFacePatch& Foam::faceZone::patch() const
{
    if (!patchPtr_)
    {
        calcFaceZonePatch();
    }
    return *patchPtr_;
}


const Foam::labelList& Foam::faceZone::frontCells() const
{
    if (!frontCellsPtr_)
    {
        calcCellLayers();
    }
    return *frontCellsPtr_;
}


const Foam::labelList& Foam::faceZone::backCells() const
{
    if (!backCellsPtr_)
    {
        calcCellLayers();
    }
    return *backCellsPtr_;
}


const Foam::labelList& Foam::faceZone::meshEdges() const
{
    if (!mePtr_)
    {
        mePtr_.reset
        (
            new labelList
            (
                this->patch().meshEdges
                (
                    zoneMesh().mesh().edges(),
                    zoneMesh().mesh().pointEdges()
                )
            )
        );
    }

    return *mePtr_;
}


void Foam::faceZone::clearGeom()
{
    patchPtr_.reset(nullptr);
    frontCellsPtr_.reset(nullptr);
    backCellsPtr_.reset(nullptr);
    mePtr_.reset(nullptr);
}


void Foam::faceZone::clearAddressing()
{
    zone::clearAddressing();
    clearGeom();
}


void Foam::faceZone::clearPrimitives()
{
    zone::clearPrimitives();
    flipMap_.clear();
}


void Foam::faceZone::resetAddressing(faceZone&& zn)
{
    // TDB: clearGeom();
    if (this == &zn)
    {
        return;  // Self-assignment is a no-op
    }

    clearAddressing();
    labelList::transfer(static_cast<labelList&>(zn));
    flipMap_.transfer(zn.flipMap_);
    zn.clearAddressing();
}


void Foam::faceZone::resetAddressing(const faceZone& zn)
{
    // TDB: clearGeom();
    if (this == &zn)
    {
        return;  // Self-assignment is a no-op
    }

    clearAddressing();
    labelList::operator=(static_cast<const labelList&>(zn));
    flipMap_ = zn.flipMap_;
}


void Foam::faceZone::resetAddressing
(
    const labelUList& addr,
    const bool flipMapValue
)
{
    clearAddressing();
    labelList::operator=(addr);
    setFlipMap(flipMapValue);
}


void Foam::faceZone::resetAddressing
(
    labelList&& addr,
    const bool flipMapValue
)
{
    clearAddressing();
    labelList::transfer(addr);
    setFlipMap(flipMapValue);
}


void Foam::faceZone::resetAddressing
(
    const labelUList& addr,
    const boolUList& flipMap
)
{
    clearAddressing();
    labelList::operator=(addr);
    flipMap_ = flipMap;
}


void Foam::faceZone::updateMesh(const mapPolyMesh& mpm)
{
    clearAddressing();

    labelList newAddressing(size());
    boolList newFlipMap(flipMap_.size());
    label nFaces = 0;

    const labelList& addr = *this;
    const boolList& flips = flipMap();
    const labelList& faceMap = mpm.reverseFaceMap();

    forAll(addr, i)
    {
        const label facei = addr[i];

        if (faceMap[facei] >= 0)
        {
            newAddressing[nFaces] = faceMap[facei];
            newFlipMap[nFaces] = flips[i];  // Keep flip map
            ++nFaces;
        }
    }

    newAddressing.resize(nFaces);
    newFlipMap.resize(nFaces);

    labelList::transfer(newAddressing);
    flipMap_.transfer(newFlipMap);
}


bool Foam::faceZone::checkDefinition(const bool report) const
{
    return zone::checkDefinition(zoneMesh().mesh().faces().size(), report);
}


bool Foam::faceZone::checkParallelSync(const bool report) const
{
    const polyMesh& mesh = zoneMesh().mesh();
    const polyBoundaryMesh& bm = mesh.boundaryMesh();

    bool hasError = false;


    // Check that zone faces are synced
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    {
        const labelList& addr = *this;
        const boolList& flips = flipMap();

        boolList neiZoneFace(mesh.nBoundaryFaces(), false);
        boolList neiZoneFlip(mesh.nBoundaryFaces(), false);

        forAll(addr, i)
        {
            const label facei = addr[i];

            if (!mesh.isInternalFace(facei))
            {
                const label bFacei = facei-mesh.nInternalFaces();
                neiZoneFace[bFacei] = true;
                neiZoneFlip[bFacei] = flips[i];
            }
        }
        boolList myZoneFace(neiZoneFace);
        boolList myZoneFlip(neiZoneFlip);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFace);
        syncTools::swapBoundaryFaceList(mesh, neiZoneFlip);

        forAll(addr, i)
        {
            const label facei = addr[i];
            const label patchi = bm.whichPatch(facei);

            if (patchi != -1 && bm[patchi].coupled())
            {
                const label bFacei = facei-mesh.nInternalFaces();

                // Check face in zone on both sides
                if (myZoneFace[bFacei] != neiZoneFace[bFacei])
                {
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << facei
                            << " on coupled patch " << bm[patchi].name()
                            << " is inconsistent with its coupled neighbour."
                            << endl;
                    }
                    else
                    {
                        // w/o report - can stop checking now
                        break;
                    }
                }
                else if (myZoneFlip[bFacei] == neiZoneFlip[bFacei])
                {
                    // Flip state should be opposite.
                    hasError = true;

                    if (report)
                    {
                        Pout<< " ***Problem with faceZone " << index()
                            << " named " << name()
                            << ". Face " << facei
                            << " on coupled patch " << bm[patchi].name()
                            << " has inconsistent flipMap across coupled faces."
                            << endl;
                    }
                    else
                    {
                        // w/o report - can stop checking now
                        break;
                    }
                }
            }
        }
    }

    return returnReduceOr(hasError);
}


void Foam::faceZone::movePoints(const pointField& pts)
{
    if (patchPtr_)
    {
        patchPtr_->movePoints(pts);
    }
}

void Foam::faceZone::write(Ostream& os) const
{
    os  << nl << name()
        << nl << static_cast<const labelList&>(*this)
        << nl << flipMap();
}


void Foam::faceZone::writeDict(Ostream& os) const
{
    os.beginBlock(name());

    os.writeEntry("type", type());
    zoneIdentifier::write(os);
    writeEntry(this->labelsName, os);
    flipMap().writeEntry("flipMap", os);

    os.endBlock();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::faceZone::operator=(const faceZone& zn)
{
    if (this == &zn)
    {
        return;  // Self-assignment is a no-op
    }

    clearAddressing();
    labelList::operator=(static_cast<const labelList&>(zn));
    flipMap_ = zn.flipMap_;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const faceZone& zn)
{
    zn.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
