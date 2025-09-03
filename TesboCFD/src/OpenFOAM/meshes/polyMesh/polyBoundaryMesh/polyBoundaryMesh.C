/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "processorPolyPatch.H"
#include "PstreamBuffers.H"
#include "lduSchedule.H"
#include "globalMeshData.H"
#include "wordRes.H"
#include "DynamicList.H"
#include "PtrListOps.H"
#include "edgeHashes.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(polyBoundaryMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::polyBoundaryMesh::hasGroupIDs() const
{
    if (groupIDsPtr_)
    {
        // Use existing cache
        return !groupIDsPtr_->empty();
    }

    const polyPatchList& patches = *this;

    for (const polyPatch& p : patches)
    {
        if (!p.inGroups().empty())
        {
            return true;
        }
    }

    return false;
}


void Foam::polyBoundaryMesh::calcGroupIDs() const
{
    if (groupIDsPtr_)
    {
        return;  // Or FatalError
    }

    groupIDsPtr_.emplace(16);
    auto& groupLookup = *groupIDsPtr_;

    const polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        for (const word& groupName : patches[patchi].inGroups())
        {
            groupLookup(groupName).push_back(patchi);
        }
    }

    // Remove groups that clash with patch names
    forAll(patches, patchi)
    {
        if (groupLookup.erase(patches[patchi].name()))
        {
            WarningInFunction
                << "Removed group '" << patches[patchi].name()
                << "' which clashes with patch " << patchi
                << " of the same name."
                << endl;
        }
    }
}


void Foam::polyBoundaryMesh::populate(PtrList<entry>&& entries)
{
    clearLocalAddressing();

    polyPatchList& patches = *this;

    patches.resize_null(entries.size());

    // Transcribe.
    // Does not handle nullptr at all (what could possibly be done?)
    forAll(patches, patchi)
    {
        patches.set
        (
            patchi,
            polyPatch::New
            (
                entries[patchi].keyword(),
                entries[patchi].dict(),
                patchi,
                *this
            )
        );
    }

    entries.clear();
}


bool Foam::polyBoundaryMesh::readIOcontents(const bool allowOptionalRead)
{
    bool updated = false;
    PtrList<entry> entries;

    if
    (
        this->isReadRequired()
     || (allowOptionalRead && this->isReadOptional() && this->headerOk())
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<polyBoundaryMesh>();

        // Read entries
        Istream& is = readStream(typeName);

        is >> entries;

        is.check(FUNCTION_NAME);
        close();
        updated = true;
    }
    // Future: support master-only and broadcast?

    if (updated)
    {
        populate(std::move(entries));
    }

    return updated;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& mesh
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(mesh)
{
    readIOcontents(false);  // allowOptionalRead = false
}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    Foam::zero
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(pm)
{}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const label size
)
:
    polyPatchList(size),
    regIOobject(io),
    mesh_(pm)
{}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    const polyPatchList& list
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(pm)
{
    if (!readIOcontents(true))  // allowOptionalRead = true
    {
        // Nothing read. Use supplied patches
        polyPatchList& patches = *this;
        patches.resize(list.size());

        forAll(patches, patchi)
        {
            patches.set(patchi, list[patchi].clone(*this));
        }
    }
}


Foam::polyBoundaryMesh::polyBoundaryMesh
(
    const IOobject& io,
    const polyMesh& pm,
    PtrList<entry>&& entries
)
:
    polyPatchList(),
    regIOobject(io),
    mesh_(pm)
{
    if (!readIOcontents(true))  // allowOptionalRead = true
    {
        populate(std::move(entries));
    }
    entries.clear();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::clear()
{
    polyPatchList::clear();
    clearAddressing();
}


void Foam::polyBoundaryMesh::clearGeom()
{
    polyPatchList& patches = *this;

    for (polyPatch& p : patches)
    {
        p.clearGeom();
    }
}


// Private until it is more generally required (and gets a better name?)
void Foam::polyBoundaryMesh::clearLocalAddressing()
{
    neighbourEdgesPtr_.reset(nullptr);
    patchIDPtr_.reset(nullptr);
    groupIDsPtr_.reset(nullptr);
}


void Foam::polyBoundaryMesh::clearAddressing()
{
    clearLocalAddressing();

    polyPatchList& patches = *this;

    for (polyPatch& p : patches)
    {
        p.clearAddressing();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMesh::calcGeometry()
{
    // Make sure messages don't interact by having unique tag
    const int oldTag = UPstream::incrMsgType();

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::buffered
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }

    // Reset tag
    UPstream::msgType(oldTag);
}


const Foam::faceList::subList Foam::polyBoundaryMesh::faces() const
{
    return faceList::subList
    (
        mesh_.faces(),
        mesh_.nBoundaryFaces(),
        mesh_.nInternalFaces()
    );
}


const Foam::labelList::subList Foam::polyBoundaryMesh::faceOwner() const
{
    return labelList::subList
    (
        mesh_.faceOwner(),
        mesh_.nBoundaryFaces(),
        mesh_.nInternalFaces()
    );
}

// Potentially useful to simplify logic elsewhere?
// const Foam::labelList::subList Foam::polyBoundaryMesh::faceNeighbour() const
// {
//     return labelList::subList();
// }


Foam::UPtrList<const Foam::labelUList>
Foam::polyBoundaryMesh::faceCells() const
{
    const polyPatchList& patches = *this;

    UPtrList<const labelUList> list(patches.size());

    forAll(patches, patchi)
    {
        list.set(patchi, &patches[patchi].faceCells());
    }

    return list;
}


const Foam::List<Foam::labelPairList>&
Foam::polyBoundaryMesh::neighbourEdges() const
{
    if (Pstream::parRun())
    {
        WarningInFunction
            << "Neighbour edge addressing not correct across parallel"
            << " boundaries." << endl;
    }

    if (!neighbourEdgesPtr_)
    {
        neighbourEdgesPtr_.emplace(size());
        auto& neighbourEdges = *neighbourEdgesPtr_;

        // Initialize.
        label nEdgePairs = 0;
        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            neighbourEdges[patchi].setSize(pp.nEdges() - pp.nInternalEdges());

            for (labelPair& edgeInfo : neighbourEdges[patchi])
            {
                edgeInfo[0] = -1;
                edgeInfo[1] = -1;
            }

            nEdgePairs += pp.nEdges() - pp.nInternalEdges();
        }

        // From mesh edge (expressed as a point pair so as not to construct
        // point addressing) to patch + relative edge index.
        EdgeMap<labelPair> pointsToEdge(nEdgePairs);

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const edgeList& edges = pp.edges();

            for
            (
                label edgei = pp.nInternalEdges();
                edgei < edges.size();
                edgei++
            )
            {
                // Edge in patch local points
                const edge& e = edges[edgei];

                // Edge in mesh points.
                edge meshEdge(pp.meshPoints()[e[0]], pp.meshPoints()[e[1]]);

                auto fnd = pointsToEdge.find(meshEdge);

                if (!fnd.good())
                {
                    // First occurrence of mesh edge. Store patch and my
                    // local index.
                    pointsToEdge.insert
                    (
                        meshEdge,
                        labelPair
                        (
                            patchi,
                            edgei - pp.nInternalEdges()
                        )
                    );
                }
                else
                {
                    // Second occurrence. Store.
                    const labelPair& edgeInfo = fnd.val();

                    neighbourEdges[patchi][edgei - pp.nInternalEdges()] =
                        edgeInfo;

                    neighbourEdges[edgeInfo[0]][edgeInfo[1]]
                         = labelPair(patchi, edgei - pp.nInternalEdges());

                    // Found all two occurrences of this edge so remove from
                    // hash to save space. Note that this will give lots of
                    // problems if the polyBoundaryMesh is multiply connected.
                    pointsToEdge.erase(meshEdge);
                }
            }
        }

        if (pointsToEdge.size())
        {
            FatalErrorInFunction
                << "Not all boundary edges of patches match up." << nl
                << "Is the outside of your mesh multiply connected?"
                << abort(FatalError);
        }

        forAll(*this, patchi)
        {
            const polyPatch& pp = operator[](patchi);

            const labelPairList& nbrEdges = neighbourEdges[patchi];

            forAll(nbrEdges, i)
            {
                const labelPair& edgeInfo = nbrEdges[i];

                if (edgeInfo[0] == -1 || edgeInfo[1] == -1)
                {
                    const label edgei = pp.nInternalEdges() + i;
                    const edge& e = pp.edges()[edgei];

                    FatalErrorInFunction
                        << "Not all boundary edges of patches match up." << nl
                        << "Edge " << edgei << " on patch " << pp.name()
                        << " end points " << pp.localPoints()[e[0]] << ' '
                        << pp.localPoints()[e[1]] << " is not matched to an"
                        << " edge on any other patch." << nl
                        << "Is the outside of your mesh multiply connected?"
                        << abort(FatalError);
                }
            }
        }
    }

    return *neighbourEdgesPtr_;
}


const Foam::labelList& Foam::polyBoundaryMesh::patchID() const
{
    if (!patchIDPtr_)
    {
        patchIDPtr_.emplace(mesh_.nBoundaryFaces());
        auto& list = *patchIDPtr_;

        const polyPatchList& patches = *this;

        forAll(patches, patchi)
        {
            SubList<label>
            (
                list,
                patches[patchi].size(),
                (patches[patchi].start() - mesh_.nInternalFaces())
            ) = patchi;
        }
    }

    return *patchIDPtr_;
}


Foam::label Foam::polyBoundaryMesh::patchID(const label meshFacei) const
{
    const label bndFacei = (meshFacei - mesh_.nInternalFaces());

    return
    (
        (bndFacei >= 0 && bndFacei < mesh_.nBoundaryFaces())
      ? this->patchID()[bndFacei]
      : -1
    );
}


Foam::labelList
Foam::polyBoundaryMesh::patchID(const labelUList& meshFaceIndices) const
{
    labelList output(meshFaceIndices.size());
    forAll(meshFaceIndices, i)
    {
        output[i] = patchID(meshFaceIndices[i]);
    }
    return output;
}


const Foam::HashTable<Foam::labelList>&
Foam::polyBoundaryMesh::groupPatchIDs() const
{
    if (!groupIDsPtr_)
    {
        calcGroupIDs();
    }

    return *groupIDsPtr_;
}


void Foam::polyBoundaryMesh::setGroup
(
    const word& groupName,
    const labelUList& patchIDs
)
{
    groupIDsPtr_.reset(nullptr);

    polyPatchList& patches = *this;

    boolList pending(patches.size(), true);

    // Add to specified patches
    for (const label patchi : patchIDs)
    {
        if (pending.test(patchi))
        {
            pending.unset(patchi);
            patches[patchi].addGroup(groupName);
        }
    }

    // Remove from other patches
    forAll(patches, patchi)
    {
        if (pending.test(patchi))
        {
            patches[patchi].removeGroup(groupName);
        }
    }
}


Foam::label Foam::polyBoundaryMesh::nNonProcessor() const
{
    const polyPatchList& patches = *this;

    label count = 0;

    for (const polyPatch& p : patches)
    {
        if (isA<processorPolyPatch>(p))
        {
            break;
        }

        ++count;
    }

    return count;
}


Foam::label Foam::polyBoundaryMesh::nProcessorPatches() const
{
    const polyPatchList& patches = *this;

    label count = 0;

    for (const polyPatch& p : patches)
    {
        if (isA<processorPolyPatch>(p))
        {
            ++count;
        }
    }

    return count;
}


Foam::wordList Foam::polyBoundaryMesh::names() const
{
    return PtrListOps::get<word>(*this, nameOp<polyPatch>());
}


Foam::wordList Foam::polyBoundaryMesh::types() const
{
    return PtrListOps::get<word>(*this, typeOp<polyPatch>());
}


Foam::wordList Foam::polyBoundaryMesh::physicalTypes() const
{
    return
        PtrListOps::get<word>
        (
            *this,
            [](const polyPatch& p) { return p.physicalType(); }
        );
}


Foam::labelList Foam::polyBoundaryMesh::patchStarts() const
{
    return
        PtrListOps::get<label>
        (
            *this,
            [](const polyPatch& p) { return p.start(); }
        );
}


Foam::labelList Foam::polyBoundaryMesh::patchSizes() const
{
    return
        PtrListOps::get<label>
        (
            *this,
            [](const polyPatch& p) { return p.size(); }
        );
}


Foam::List<Foam::labelRange> Foam::polyBoundaryMesh::patchRanges() const
{
    return
        PtrListOps::get<labelRange>
        (
            *this,
            [](const polyPatch& p) { return p.range(); }
        );
}


Foam::wordList Foam::polyBoundaryMesh::groupNames() const
{
    return this->groupPatchIDs().sortedToc();
}


Foam::label Foam::polyBoundaryMesh::start() const noexcept
{
    return mesh_.nInternalFaces();
}


Foam::label Foam::polyBoundaryMesh::nFaces() const noexcept
{
    return mesh_.nBoundaryFaces();
}


Foam::labelRange Foam::polyBoundaryMesh::range() const noexcept
{
    return labelRange(mesh_.nInternalFaces(), mesh_.nBoundaryFaces());
}


Foam::labelRange Foam::polyBoundaryMesh::range(const label patchi) const
{
    if (patchi < 0)
    {
        return labelRange(mesh_.nInternalFaces(), 0);
    }

    // Will fail if patchi >= size()
    return (*this)[patchi].range();
}


Foam::labelList Foam::polyBoundaryMesh::indices
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    if (matcher.empty())
    {
        return labelList();
    }

    // Only check groups if requested and they exist
    const bool checkGroups = (useGroups && this->hasGroupIDs());

    labelHashSet ids(0);

    if (matcher.isPattern())
    {
        if (checkGroups)
        {
            const auto& groupLookup = groupPatchIDs();
            forAllConstIters(groupLookup, iter)
            {
                if (matcher(iter.key()))
                {
                    // Add patch ids associated with the group
                    ids.insert(iter.val());
                }
            }
        }

        if (ids.empty())
        {
            return PtrListOps::findMatching(*this, matcher);
        }
        else
        {
            ids.insert(PtrListOps::findMatching(*this, matcher));
        }
    }
    else
    {
        // Literal string.
        // Special version of above for reduced memory footprint.

        const label patchId = PtrListOps::firstMatching(*this, matcher);

        if (patchId >= 0)
        {
            return labelList(one{}, patchId);
        }
        else if (checkGroups)
        {
            const auto iter = groupPatchIDs().cfind(matcher);

            if (iter.good())
            {
                // Hash ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    return ids.sortedToc();
}


Foam::labelList Foam::polyBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    if (matcher.empty())
    {
        return labelList();
    }
    else if (matcher.size() == 1)
    {
        return this->indices(matcher.front(), useGroups);
    }

    labelHashSet ids(0);

    // Only check groups if requested and they exist
    if (useGroups && this->hasGroupIDs())
    {
        ids.reserve(this->size());

        const auto& groupLookup = groupPatchIDs();
        forAllConstIters(groupLookup, iter)
        {
            if (matcher(iter.key()))
            {
                // Add patch ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    if (ids.empty())
    {
        return PtrListOps::findMatching(*this, matcher);
    }
    else
    {
        ids.insert(PtrListOps::findMatching(*this, matcher));
    }

    return ids.sortedToc();
}


Foam::labelList Foam::polyBoundaryMesh::indices
(
    const wordRes& select,
    const wordRes& ignore,
    const bool useGroups
) const
{
    if (ignore.empty())
    {
        return this->indices(select, useGroups);
    }

    const wordRes::filter matcher(select, ignore);

    labelHashSet ids(0);

    // Only check groups if requested and they exist
    if (useGroups && this->hasGroupIDs())
    {
        ids.reserve(this->size());

        const auto& groupLookup = groupPatchIDs();
        forAllConstIters(groupLookup, iter)
        {
            if (matcher(iter.key()))
            {
                // Add patch ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    if (ids.empty())
    {
        return PtrListOps::findMatching(*this, matcher);
    }
    else
    {
        ids.insert(PtrListOps::findMatching(*this, matcher));
    }

    return ids.sortedToc();
}


Foam::label Foam::polyBoundaryMesh::findIndex(const wordRe& key) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


Foam::label Foam::polyBoundaryMesh::findPatchID
(
    const word& patchName,
    bool allowNotFound
) const
{
    if (patchName.empty())
    {
        return -1;
    }

    const label patchId = PtrListOps::firstMatching(*this, patchName);

    if (patchId >= 0)
    {
        return patchId;
    }

    if (!allowNotFound)
    {
        FatalErrorInFunction
            << "Patch '" << patchName << "' not found. "
            << "Available patch names";

        if (polyMesh::defaultRegion != mesh_.name())
        {
            FatalError
                << " in region '" << mesh_.name() << "'";
        }

        FatalError
            << " include: " << names() << endl
            << exit(FatalError);
    }

    // Patch not found
    if (debug)
    {
        Pout<< "label polyBoundaryMesh::findPatchID(const word&) const"
            << "Patch named " << patchName << " not found.  "
            << "Available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


Foam::labelPair
Foam::polyBoundaryMesh::whichPatchFace(const label meshFacei) const
{
    if (meshFacei < mesh().nInternalFaces())
    {
        // Internal face: return (-1, meshFace)
        return labelPair(-1, meshFacei);
    }
    else if (meshFacei >= mesh().nFaces())
    {
        // Bounds error: abort
        FatalErrorInFunction
            << "Face " << meshFacei
            << " out of bounds. Number of geometric faces " << mesh().nFaces()
            << abort(FatalError);

        return labelPair(-1, meshFacei);
    }


    const polyPatchList& patches = *this;

    // Do we have cached patch info?
    if (patchIDPtr_)
    {
        const label patchi =
            this->patchID()[meshFacei - mesh().nInternalFaces()];

        // (patch, local face index)
        return labelPair(patchi, meshFacei - patches[patchi].start());
    }


    // Patches are ordered, use binary search
    // Find out which patch index and local patch face the specified
    // mesh face belongs to by comparing label with patch start labels.

    const label patchi =
        findLower
        (
            patches,
            meshFacei,
            0,
            // Must include the start in the comparison
            [](const polyPatch& p, label val) { return (p.start() <= val); }
        );

    if (patchi < 0 || !patches[patchi].range().contains(meshFacei))
    {
        // If not in any of above, it is trouble!
        FatalErrorInFunction
            << "Face " << meshFacei << " not found in any of the patches "
            << flatOutput(names()) << nl
            << "The patches appear to be inconsistent with the mesh :"
            << " internalFaces:" << mesh().nInternalFaces()
            << " total number of faces:" << mesh().nFaces()
            << abort(FatalError);

        return labelPair(-1, meshFacei);
    }

    // (patch, local face index)
    return labelPair(patchi, meshFacei - patches[patchi].start());
}


Foam::labelPairList
Foam::polyBoundaryMesh::whichPatchFace(const labelUList& meshFaceIndices) const
{
    labelPairList output(meshFaceIndices.size());
    forAll(meshFaceIndices, i)
    {
        output[i] = whichPatchFace(meshFaceIndices[i]);
    }
    return output;
}


Foam::labelHashSet Foam::polyBoundaryMesh::patchSet
(
    const UList<wordRe>& select,
    const bool warnNotFound,
    const bool useGroups
) const
{
    labelHashSet ids(0);
    if (select.empty())
    {
        return ids;
    }

    const polyPatchList& patches = *this;

    const label len = patches.size();

    ids.reserve(len);

    // Only check groups if requested and they exist
    const bool checkGroups = (useGroups && this->hasGroupIDs());

    for (const wordRe& matcher : select)
    {
        bool missed = true;

        for (label i = 0; i < len; ++i)
        {
            if (matcher(patches[i].name()))
            {
                ids.insert(i);
                missed = false;
            }
        }

        if (missed && checkGroups)
        {
            // Check group names
            if (matcher.isPattern())
            {
                forAllConstIters(groupPatchIDs(), iter)
                {
                    if (matcher.match(iter.key()))
                    {
                        // Hash ids associated with the group
                        ids.insert(iter.val());
                        missed = false;
                    }
                }
            }
            else
            {
                const auto iter = groupPatchIDs().cfind(matcher);

                if (iter.good())
                {
                    // Hash ids associated with the group
                    ids.insert(iter.val());
                    missed = false;
                }
            }
        }

        if (missed && warnNotFound)
        {
            if (checkGroups)
            {
                WarningInFunction
                    << "Cannot find any patch or group names matching "
                    << matcher << endl;
            }
            else
            {
                WarningInFunction
                    << "Cannot find any patch names matching "
                    << matcher << endl;
            }
        }
    }

    return ids;
}


void Foam::polyBoundaryMesh::matchGroups
(
    const labelUList& patchIDs,
    wordList& groups,
    labelHashSet& nonGroupPatches
) const
{
    // Current matched groups
    DynamicList<word> matchedGroups(1);

    // Current set of unmatched patches
    nonGroupPatches = labelHashSet(patchIDs);

    const HashTable<labelList>& groupLookup = this->groupPatchIDs();
    forAllConstIters(groupLookup, iter)
    {
        // Store currently unmatched patches so we can restore
        labelHashSet oldNonGroupPatches(nonGroupPatches);

        // Match by deleting patches in group from the current set and seeing
        // if all have been deleted.
        labelHashSet groupPatchSet(iter.val());

        label nMatch = nonGroupPatches.erase(groupPatchSet);

        if (nMatch == groupPatchSet.size())
        {
            matchedGroups.push_back(iter.key());
        }
        else if (nMatch != 0)
        {
            // No full match. Undo.
            nonGroupPatches.transfer(oldNonGroupPatches);
        }
    }

    groups.transfer(matchedGroups);
}


bool Foam::polyBoundaryMesh::checkParallelSync(const bool report) const
{
    if (!Pstream::parRun())
    {
        return false;
    }

    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    // Collect non-proc patches and check proc patches are last.
    wordList localNames(bm.size());
    wordList localTypes(bm.size());

    label nonProci = 0;

    forAll(bm, patchi)
    {
        if (!isA<processorPolyPatch>(bm[patchi]))
        {
            if (nonProci != patchi)
            {
                // A processor patch in between normal patches!
                hasError = true;

                if (debug || report)
                {
                    Pout<< " ***Problem with boundary patch " << patchi
                        << " name:" << bm[patchi].name()
                        << " type:" <<  bm[patchi].type()
                        << " - seems to be preceeded by processor patches."
                        << " This is usually a problem." << endl;
                }
            }
            else
            {
                localNames[nonProci] = bm[patchi].name();
                localTypes[nonProci] = bm[patchi].type();
                ++nonProci;
            }
        }
    }
    localNames.resize(nonProci);
    localTypes.resize(nonProci);

    // Check and report error(s) on master
    // - don't need indexing on master itself

    const globalIndex procAddr(globalIndex::gatherNonLocal{}, nonProci);

    const wordList allNames(procAddr.gather(localNames));
    const wordList allTypes(procAddr.gather(localTypes));

    // Automatically restricted to master
    for (const int proci : procAddr.subProcs())
    {
        const auto procNames(allNames.slice(procAddr.range(proci)));
        const auto procTypes(allTypes.slice(procAddr.range(proci)));

        if (procNames != localNames || procTypes != localTypes)
        {
            hasError = true;

            if (debug || report)
            {
                Info<< " ***Inconsistent patches across processors, "
                       "processor0 has patch names:" << localNames
                    << " patch types:" << localTypes
                    << " processor" << proci
                    << " has patch names:" << procNames
                    << " patch types:" << procTypes
                    << endl;
            }
        }
    }

    // Reduce (not broadcast) to respect local out-of-order errors (first loop)
    return returnReduceOr(hasError);
}


bool Foam::polyBoundaryMesh::checkDefinition(const bool report) const
{
    label nextPatchStart = mesh().nInternalFaces();
    const polyBoundaryMesh& bm = *this;

    bool hasError = false;

    wordHashSet patchNames(2*this->size());

    forAll(bm, patchi)
    {
        if (bm[patchi].start() != nextPatchStart && !hasError)
        {
            hasError = true;

            Info<< " ****Problem with boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << ". The patch should start on face no " << nextPatchStart
                << " and the patch specifies " << bm[patchi].start()
                << "." << endl
                << "Possibly consecutive patches have this same problem."
                << " Suppressing future warnings." << endl;
        }

        if (!patchNames.insert(bm[patchi].name()) && !hasError)
        {
            hasError = true;

            Info<< " ****Duplicate boundary patch " << patchi
                << " named " << bm[patchi].name()
                << " of type " <<  bm[patchi].type()
                << "." << endl
                << "Suppressing future warnings." << endl;
        }

        nextPatchStart += bm[patchi].size();
    }

    Pstream::reduceOr(hasError);

    if (debug || report)
    {
        if (hasError)
        {
            Pout<< " ***Boundary definition is in error." << endl;
        }
        else
        {
            Info<< "    Boundary definition OK." << endl;
        }
    }

    return hasError;
}


void Foam::polyBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::buffered
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::polyBoundaryMesh::updateMesh()
{
    neighbourEdgesPtr_.reset(nullptr);
    patchIDPtr_.reset(nullptr);
    groupIDsPtr_.reset(nullptr);

    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::buffered
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initUpdateMesh(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).updateMesh(pBufs);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initUpdateMesh(pBufs);
            }
            else
            {
                operator[](patchi).updateMesh(pBufs);
            }
        }
    }
}


void Foam::polyBoundaryMesh::reorder
(
    const labelUList& oldToNew,
    const bool validBoundary
)
{
    // Change order of patches
    polyPatchList::reorder(oldToNew);

    // Adapt indices
    polyPatchList& patches = *this;

    forAll(patches, patchi)
    {
        patches[patchi].index() = patchi;
    }

    // Clear group-to-patch addressing. Note: could re-calculate
    groupIDsPtr_.reset(nullptr);

    if (validBoundary)
    {
        updateMesh();
    }
}


void Foam::polyBoundaryMesh::writeEntry(Ostream& os) const
{
    const polyPatchList& entries = *this;

    os  << entries.size();

    if (entries.empty())
    {
        // 0-sized : can write with less vertical space
        os  << token::BEGIN_LIST << token::END_LIST;
    }
    else
    {
        os  << nl << token::BEGIN_LIST << incrIndent << nl;

        for (const auto& pp : entries)
        {
            os.beginBlock(pp.name());
            os  << pp;
            os.endBlock();
        }
        os  << decrIndent << token::END_LIST;
    }
    os.check(FUNCTION_NAME);
}


void Foam::polyBoundaryMesh::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    const polyPatchList& entries = *this;

    if (!keyword.empty())
    {
        os.write(keyword);
        os << (entries.empty() ? token::SPACE : token::NL);
    }

    writeEntry(os);

    if (!keyword.empty()) os.endEntry();
}


bool Foam::polyBoundaryMesh::writeData(Ostream& os) const
{
    writeEntry(os);
    return os.good();
}


bool Foam::polyBoundaryMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    streamOpt.compression(IOstreamOption::UNCOMPRESSED);
    return regIOobject::writeObject(streamOpt, writeOnProc);
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

const Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
) const
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


Foam::polyPatch& Foam::polyBoundaryMesh::operator[]
(
    const word& patchName
)
{
    const label patchi = findPatchID(patchName);

    if (patchi < 0)
    {
        FatalErrorInFunction
            << "Patch named " << patchName << " not found." << nl
            << "Available patch names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](patchi);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyBoundaryMesh& pbm)
{
    pbm.writeData(os);
    return os;
}


// ************************************************************************* //
