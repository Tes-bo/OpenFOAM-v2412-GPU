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

#include "pointBoundaryMesh.H"
#include "polyBoundaryMesh.H"
#include "facePointPatch.H"
#include "pointMesh.H"
#include "PstreamBuffers.H"
#include "lduSchedule.H"
#include "globalMeshData.H"
#include "processorPointPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointBoundaryMesh, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::pointBoundaryMesh::hasGroupIDs() const
{
    if (groupIDsPtr_)
    {
        // Use existing cache
        return !groupIDsPtr_->empty();
    }

    const auto& patches = *this;

    for (const auto& p : patches)
    {
        if (!p.inGroups().empty())
        {
            return true;
        }
    }

    return false;
}


void Foam::pointBoundaryMesh::calcGroupIDs() const
{
    if (groupIDsPtr_)
    {
        return;  // Or FatalError
    }

    groupIDsPtr_.emplace(16);
    auto& groupLookup = *groupIDsPtr_;

    const auto& patches = *this;

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


void Foam::pointBoundaryMesh::addPatches(const polyBoundaryMesh& pbm)
{
    // Set boundary patches
    pointPatchList& patches = *this;

    patches.resize_null(pbm.size());

    forAll(pbm, patchi)
    {
        // NB: needs ptr() to get *pointPatch instead of *facePointPatch
        patches.set(patchi, facePointPatch::New(pbm[patchi], *this).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointBoundaryMesh::pointBoundaryMesh
(
    const pointMesh& m,
    const polyBoundaryMesh& pbm
)
:
    pointPatchList(),
    regIOobject
    (
        IOobject
        (
            "boundary",
            //m.thisDb().time().findInstance(m.meshDir(), "boundary"),
            pbm.mesh().facesInstance(),
            polyMesh::meshSubDir/pointMesh::meshSubDir,
            m.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false                   // Avoid conflict with polyMesh/boundary
        )
    ),
    mesh_(m)
{
    addPatches(pbm);

    if (debug)
    {
        pointPatchList& Patches = *this;

        Pout<< "pointBoundaryMesh::pointBoundaryMesh"
            << "(const pointMesh&, const polyBoundaryMesh&): "
            << "constructed pointBoundaryMesh:" << endl;
        Pout<< incrIndent;
        for (const auto& pp : Patches)
        {
            Pout<< indent
                << "index:" << pp.index() << " patch:" << pp.name()
                << " type:" << pp.type() << endl;
        }
        Pout<< decrIndent;
    }
}


Foam::pointBoundaryMesh::pointBoundaryMesh
(
    const IOobject& io,
    const pointMesh& m,
    const polyBoundaryMesh& pbm
)
:
    pointPatchList(),
    regIOobject
    (
        IOobject
        (
            "boundary",
            io.instance(),
            polyMesh::meshSubDir/pointMesh::meshSubDir,
            io.db(),
            io.readOpt(),
            io.writeOpt(),
            false   //io.registerObject()     // or always set to false?
        )
    ),
    mesh_(m)
{
    pointPatchList& Patches = *this;

    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        if (readOpt() == IOobject::MUST_READ_IF_MODIFIED)
        {
            WarningInFunction
                << "Specified IOobject::MUST_READ_IF_MODIFIED but class"
                << " does not support automatic rereading."
                << endl;
        }

        if (debug)
        {
            Pout<< "pointBoundaryMesh::pointBoundaryMesh"
                << "(const IOobject&, const pointMesh&,"
                << " const polyBoundaryMesh&): "
                << "Constructing from boundary file " << objectRelPath()
                << endl;
        }

        // Read pointPatchList
        Istream& is = readStream(typeName);

        PtrList<entry> patchEntries(is);
        Patches.setSize(patchEntries.size());

        forAll(Patches, patchi)
        {
            // Try construct-from-dictionary first
            const word& name = patchEntries[patchi].keyword();

            autoPtr<pointPatch> pPtr
            (
                pointPatch::New
                (
                    name,
                    patchEntries[patchi].dict(),
                    patchi,
                    *this
                )
            );

            if (!pPtr)
            {
                const label polyPatchi = pbm.findPatchID(name, false);
                // Try as facePointPatch from polyPatch
                pPtr = facePointPatch::New(pbm[polyPatchi], *this);
                pPtr->index() = patchi;
            }

            Patches.set(patchi, pPtr);
        }

        // Check state of IOstream
        is.check
        (
            "pointBoundaryMesh::pointBoundaryMesh"
            "(const IOobject&, const pointMesh&,"
            " const polyBoundaryMesh&)"
        );

        close();
    }
    else
    {
        if (debug)
        {
            Pout<< "pointBoundaryMesh::pointBoundaryMesh"
                << "(const IOobject&, const pointMesh&,"
                << " const polyBoundaryMesh&): "
                << "Constructing from polyBoundaryMesh only"
                << endl;
        }

        addPatches(pbm);
    }


    if (debug)
    {
        Pout<< "pointBoundaryMesh::pointBoundaryMesh"
            << "(const IOobject&, const pointMesh&, const polyBoundaryMesh&): "
            << "constructed pointBoundaryMesh:" << endl;
        Pout<< incrIndent;
        for (const auto& pp : Patches)
        {
            Pout<< indent
                << "index:" << pp.index() << " patch:" << pp.name()
                << " type:" << pp.type() << endl;
        }
        Pout<< decrIndent;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::pointBoundaryMesh::nNonProcessor() const
{
    const pointPatchList& patches = *this;

    label count = 0;

    for (const auto& p : patches)
    {
        if (isA<processorPointPatch>(p))
        {
            break;
        }

        ++count;
    }

    return count;
}


Foam::label Foam::pointBoundaryMesh::nProcessorPatches() const
{
    const pointPatchList& patches = *this;

    label count = 0;

    for (const auto& p : patches)
    {
        if (isA<processorPointPatch>(p))
        {
            ++count;
        }
    }

    return count;
}


Foam::wordList Foam::pointBoundaryMesh::names() const
{
    return PtrListOps::get<word>(*this, nameOp<pointPatch>());
}


Foam::wordList Foam::pointBoundaryMesh::types() const
{
    return PtrListOps::get<word>(*this, typeOp<pointPatch>());
}


Foam::wordList Foam::pointBoundaryMesh::physicalTypes() const
{
    return
        PtrListOps::get<word>
        (
            *this,
            [](const pointPatch& p) { return p.physicalType(); }
        );
}


const Foam::HashTable<Foam::labelList>&
Foam::pointBoundaryMesh::groupPatchIDs() const
{
    if (!groupIDsPtr_)
    {
        calcGroupIDs();
    }

    return *groupIDsPtr_;
}


Foam::labelList Foam::pointBoundaryMesh::indices
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


Foam::labelList Foam::pointBoundaryMesh::indices
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


Foam::labelList Foam::pointBoundaryMesh::indices
(
    const wordRes& select,
    const wordRes& ignore,
    const bool useGroups
) const
{
    //return mesh().boundaryMesh().indices(select, ignore, useGroups);
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


Foam::label Foam::pointBoundaryMesh::findPatchID
(
    const word& patchName,
    bool allowNotFound
) const
{
    //return mesh().boundaryMesh().findPatchID(patchName);
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
        Pout<< "label pointBoundaryMesh::findPatchID(const word&) const"
            << " Patch named " << patchName << " not found.  "
            << "Available patch names: " << names() << endl;
    }

    // Not found, return -1
    return -1;
}


void Foam::pointBoundaryMesh::calcGeometry()
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
}


void Foam::pointBoundaryMesh::movePoints(const pointField& p)
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


void Foam::pointBoundaryMesh::updateMesh()
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


void Foam::pointBoundaryMesh::reorder
(
    const labelUList& oldToNew,
    const bool validBoundary
)
{
    // Change order of patches
    pointPatchList::reorder(oldToNew);

    // Adapt indices
    pointPatchList& patches = *this;

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

    if (debug)
    {
        pointPatchList& Patches = *this;

        Pout<< "pointBoundaryMesh::reorder"
            << "(const labelUList&, const bool): "
            << "reordered pointBoundaryMesh:" << endl;
        Pout<< incrIndent;
        for (const auto& pp : Patches)
        {
            Pout<< indent
                << "index:" << pp.index() << " patch:" << pp.name()
                << " type:" << pp.type() << endl;
        }
        Pout<< decrIndent;
    }
}


bool Foam::pointBoundaryMesh::writeData(Ostream& os) const
{
    const pointPatchList& patches = *this;

    os  << patches.size() << nl << token::BEGIN_LIST << incrIndent << nl;

    forAll(patches, patchi)
    {
        os  << indent << patches[patchi].name() << nl
            << indent << token::BEGIN_BLOCK << nl
            << incrIndent << patches[patchi] << decrIndent
            << indent << token::END_BLOCK << endl;
    }

    os  << decrIndent << token::END_LIST;

    // Check state of IOstream
    os.check("pointBoundaryMesh::writeData(Ostream& os) const");

    return os.good();
}


//bool Foam::pointBoundaryMesh::writeObject
//(
//    IOstreamOption
//) const
//{
//    return regIOobject::writeObject(fmt, ver, IOstream::UNCOMPRESSED);
//}


// ************************************************************************* //
