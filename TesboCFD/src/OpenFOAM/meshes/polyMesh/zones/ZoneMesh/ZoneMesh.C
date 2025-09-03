/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2024 OpenCFD Ltd.
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

#include "ZoneMesh.H"
#include "entry.H"
#include "DynamicList.H"
#include "Pstream.H"
#include "PtrListOps.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    template<class ZoneType, class MeshType>
    int Foam::ZoneMesh<ZoneType, MeshType>::disallowGenericZones
    (
        debug::debugSwitch("disallowGenericZones", 0)
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::totalSize() const
{
    // Count number of objects in all zones
    const PtrList<ZoneType>& zones = *this;

    label total = 0;
    for (const ZoneType& zn : zones)
    {
        total += zn.size();
    }

    return total;
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::calcZoneMap() const
{
    if (zoneMapPtr_)
    {
        FatalErrorInFunction
            << "zone map already calculated"
            << abort(FatalError);
    }
    else
    {
        zoneMapPtr_.reset(new Map<label>());
        auto& map = *zoneMapPtr_;

        // Fill in objects of all zones into the map.
        // The key is the global object index, value is the zone index

        map.reserve(this->totalSize());

        const PtrList<ZoneType>& zones = *this;

        label zonei = 0;

        for (const ZoneType& zn : zones)
        {
            for (const label id : static_cast<const labelList&>(zn))
            {
                //map.insert(id, zonei);
                const auto fnd = map.cfind(id);
                if (!fnd)
                {
                    map.insert(id, zonei);
                }
                else if (fnd.val() != zonei)
                {
                    // Multiple zones for same id

                    if (!additionalMapPtr_)
                    {
                        // First occurrence
                        label maxIndex = -1;
                        for (const ZoneType& zn : zones)
                        {
                            for
                            (
                                const label id
                              : static_cast<const labelList&>(zn)
                            )
                            {
                                maxIndex = Foam::max(maxIndex, id);
                            }
                        }
                        additionalMapPtr_.reset(new labelListList(maxIndex+1));
                    }
                    auto& additionalMap = *additionalMapPtr_;
                    additionalMap[id].push_uniq(zonei);
                }
            }

            ++zonei;
        }

        // Sort such that map contains lowest zoneID
        if (additionalMapPtr_)
        {
            auto& additionalMap = *additionalMapPtr_;
            forAll(additionalMap, id)
            {
                labelList& zones = additionalMap[id];

                if (zones.size())
                {
                    Foam::stableSort(zones);
                    const label zonei = map[id];
                    const label index = findLower(zones, zonei);
                    if (index == -1)
                    {
                        // Already first
                    }
                    else
                    {
                        map.set(id, zones[0]);
                        for (label i = 0; i < zones.size() && i <= index; i++)
                        {
                            zones[i] = zones[i+1];
                        }
                        zones[index+1] = zonei;
                    }
                }
            }
        }
    }
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::hasGroupIDs() const
{
    if (groupIDsPtr_)
    {
        // Use existing cache
        return !groupIDsPtr_->empty();
    }

    const PtrList<ZoneType>& zones = *this;

    for (const ZoneType& zn : zones)
    {
        if (!zn.inGroups().empty())
        {
            return true;
        }
    }

    return false;
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::calcGroupIDs() const
{
    if (groupIDsPtr_)
    {
        return;  // Or FatalError
    }

    groupIDsPtr_.reset(new HashTable<labelList>(16));
    auto& groupLookup = *groupIDsPtr_;

    const PtrList<ZoneType>& zones = *this;

    forAll(zones, zonei)
    {
        for (const word& groupName : zones[zonei].inGroups())
        {
            groupLookup(groupName).push_back(zonei);
        }
    }

    // Remove groups that clash with zone names
    forAll(zones, zonei)
    {
        if (groupLookup.erase(zones[zonei].name()))
        {
            WarningInFunction
                << "Removed group '" << zones[zonei].name()
                << "' which clashes with zone " << zonei
                << " of the same name."
                << endl;
        }
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::populate
(
    PtrList<entry>&& entries
)
{
    clearLocalAddressing();

    PtrList<ZoneType>& zones = *this;

    zones.resize_null(entries.size());

    // Transcribe
    // Does not handle nullptr at all
    forAll(zones, zonei)
    {
        // Possible handling for nullptr:
        // zones.emplace_set
        // (
        //     zonei,
        //     "missing_" + ::Foam::name(zonei), zonei, *this
        // );

        zones.set
        (
            zonei,
            ZoneType::New
            (
                entries[zonei].keyword(),
                entries[zonei].dict(),
                zonei,
                *this
            )
        );
    }

    entries.clear();
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::readIOcontents
(
    const bool allowOptionalRead
)
{
    bool updated = false;
    PtrList<entry> entries;

    if
    (
        isReadRequired()
     || (allowOptionalRead && isReadOptional() && headerOk())
    )
    {
        // Warn for MUST_READ_IF_MODIFIED
        warnNoRereading<ZoneMesh<ZoneType, MeshType>>();

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

template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh)
{
    // Note: this is inconsistent with polyBoundaryMesh
    // which does not permit optional reading
    readIOcontents(true);  // allowOptionalRead = true
}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    Foam::zero
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh)
{}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    const label size
)
:
    PtrList<ZoneType>(size),
    regIOobject(io),
    mesh_(mesh)
{
    // Note: this is inconsistent with polyBoundaryMesh
    // which does not read all
    readIOcontents(true);  // allowOptionalRead = true
}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    const PtrList<ZoneType>& list
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh)
{
    if (!readIOcontents(true))  // allowOptionalRead = true
    {
        // Nothing read. Use supplied zones
        PtrList<ZoneType>& zones = *this;
        zones.resize(list.size());

        forAll(zones, zonei)
        {
            zones.set(zonei, list[zonei].clone(*this));
        }
    }
}


template<class ZoneType, class MeshType>
Foam::ZoneMesh<ZoneType, MeshType>::ZoneMesh
(
    const IOobject& io,
    const MeshType& mesh,
    PtrList<entry>&& entries
)
:
    PtrList<ZoneType>(),
    regIOobject(io),
    mesh_(mesh)
{
    if (!readIOcontents(true))  // allowOptionalRead = true
    {
        populate(std::move(entries));
    }
    entries.clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const Foam::Map<Foam::label>&
Foam::ZoneMesh<ZoneType, MeshType>::zoneMap() const
{
    if (!zoneMapPtr_)
    {
        calcZoneMap();
    }

    return *zoneMapPtr_;
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::whichZone
(
    const label objectIndex
) const
{
    return zoneMap().lookup(objectIndex, -1);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::whichZones
(
    const label objectIndex,
    DynamicList<label>& zones
) const
{
    zones.clear();
    const auto fnd = zoneMap().cfind(objectIndex);
    if (fnd)
    {
        // Add main element
        zones.push_back(fnd.val());
        if (additionalMapPtr_)
        {
            const auto& additionalMap = *additionalMapPtr_;
            if (objectIndex < additionalMap.size())
            {
                for (const label zonei : additionalMap[objectIndex])
                {
                    zones.push_uniq(zonei);
                }
            }
        }
    }
    return zones.size();
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::types() const
{
    return PtrListOps::get<word>(*this, typeOp<ZoneType>());
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names() const
{
    return PtrListOps::get<word>(*this, nameOp<ZoneType>());
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::groupNames() const
{
    return this->groupZoneIDs().sortedToc();
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names
(
    const wordRe& matcher
) const
{
    return PtrListOps::names(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::names
(
    const wordRes& matcher
)
const
{
    return PtrListOps::names(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames() const
{
    wordList sorted(this->names());
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames
(
    const wordRe& matcher
) const
{
    wordList sorted(this->names(matcher));
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::wordList Foam::ZoneMesh<ZoneType, MeshType>::sortedNames
(
    const wordRes& matcher
)
const
{
    wordList sorted(this->names(matcher));
    Foam::sort(sorted);

    return sorted;
}


template<class ZoneType, class MeshType>
Foam::labelList Foam::ZoneMesh<ZoneType, MeshType>::indices
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

    if (checkGroups)
    {
        ids.reserve(this->size());
    }

    if (matcher.isPattern())
    {
        if (checkGroups)
        {
            const auto& groupLookup = groupZoneIDs();
            forAllConstIters(groupLookup, iter)
            {
                if (matcher(iter.key()))
                {
                    // Hash ids associated with the group
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
        // Special version of above for reduced memory footprint

        const label zoneId = PtrListOps::firstMatching(*this, matcher);

        if (zoneId >= 0)
        {
            return labelList(one{}, zoneId);
        }
        else if (checkGroups)
        {
            const auto iter = groupZoneIDs().cfind(matcher);

            if (iter.good())
            {
                // Hash ids associated with the group
                ids.insert(iter.val());
            }
        }
    }

    return ids.sortedToc();
}


template<class ZoneType, class MeshType>
Foam::labelList Foam::ZoneMesh<ZoneType, MeshType>::indices
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

        const auto& groupLookup = groupZoneIDs();
        forAllConstIters(groupLookup, iter)
        {
            if (matcher(iter.key()))
            {
                // Hash the ids associated with the group
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


template<class ZoneType, class MeshType>
Foam::labelList Foam::ZoneMesh<ZoneType, MeshType>::indices
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

        const auto& groupLookup = groupZoneIDs();
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


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findIndex
(
    const wordRe& key
) const
{
    if (key.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, key);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findIndex
(
    const wordRes& matcher
) const
{
    if (matcher.empty())
    {
        return -1;
    }
    return PtrListOps::firstMatching(*this, matcher);
}


template<class ZoneType, class MeshType>
Foam::label Foam::ZoneMesh<ZoneType, MeshType>::findZoneID
(
    const word& zoneName
) const
{
    if (zoneName.empty())
    {
        return -1;
    }

    label zoneId = PtrListOps::firstMatching(*this, zoneName);

    if (zoneId < 0)
    {
        DebugInFunction
            << "Zone named " << zoneName << " not found.  "
            << "List of available zone names: " << names() << endl;

        // Used for -dry-run, for example
        if (disallowGenericZones != 0)
        {
            Info<< "Creating dummy zone " << zoneName << endl;
            auto& zm = const_cast<ZoneMesh<ZoneType, MeshType>&>(*this);
            zm.emplace_back(zoneName, zm.size(), zm);
        }
    }

    return zoneId;
}


template<class ZoneType, class MeshType>
const ZoneType* Foam::ZoneMesh<ZoneType, MeshType>::cfindZone
(
    const word& zoneName
) const
{
    if (zoneName.empty())
    {
        return nullptr;
    }

    const PtrList<ZoneType>& zones = *this;

    for (auto iter = zones.begin(); iter != zones.end(); ++iter)
    {
        const ZoneType* ptr = iter.get();

        if (ptr && zoneName == ptr->name())
        {
            return ptr;
        }
    }

    // Used for -dry-run, for example
    if (disallowGenericZones != 0)
    {
        Info<< "Creating dummy zone " << zoneName << endl;
        auto& zm = const_cast<ZoneMesh<ZoneType, MeshType>&>(*this);
        zm.emplace_back(zoneName, zm.size(), zm);
    }

    return nullptr;
}


template<class ZoneType, class MeshType>
ZoneType* Foam::ZoneMesh<ZoneType, MeshType>::findZone
(
    const word& zoneName
)
{
    return const_cast<ZoneType*>(this->cfindZone(zoneName));
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const labelUList& zoneIds
) const
{
    bitSet bitset;

    for (const label zonei : zoneIds)
    {
        #ifdef FULLDEBUG
        if (zonei < 0 || zonei >= this->size())
        {
            FatalErrorInFunction
                << ZoneType::typeName << " "
                << zonei << " out of range [0," << this->size() << ")"
                << abort(FatalError);
        }
        #endif

        bitset.set
        (
            static_cast<const labelList&>(this->operator[](zonei))
        );
    }

    return bitset;
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    // matcher.empty() is handled by indices()
    return this->selection(this->indices(matcher, useGroups));
}


template<class ZoneType, class MeshType>
Foam::bitSet Foam::ZoneMesh<ZoneType, MeshType>::selection
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    // matcher.empty() is handled by indices()
    return this->selection(this->indices(matcher, useGroups));
}


template<class ZoneType, class MeshType>
const Foam::HashTable<Foam::labelList>&
Foam::ZoneMesh<ZoneType, MeshType>::groupZoneIDs() const
{
    if (!groupIDsPtr_)
    {
        calcGroupIDs();
    }

    return *groupIDsPtr_;
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::setGroup
(
    const word& groupName,
    const labelUList& zoneIDs
)
{
    groupIDsPtr_.reset(nullptr);

    PtrList<ZoneType>& zones = *this;

    boolList pending(zones.size(), true);

    // Add to specified zones
    for (const label zonei : zoneIDs)
    {
        if (pending.test(zonei))
        {
            pending.unset(zonei);
            zones[zonei].addGroup(groupName);
        }
    }

    // Remove from other zones
    forAll(zones, zonei)
    {
        if (pending.test(zonei))
        {
            zones[zonei].removeGroup(groupName);
        }
    }
}


// Private until it is more generally required (and gets a better name?)
template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clearLocalAddressing()
{
    zoneMapPtr_.reset(nullptr);
    additionalMapPtr_.reset(nullptr);
    groupIDsPtr_.reset(nullptr);
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clearAddressing()
{
    clearLocalAddressing();

    PtrList<ZoneType>& zones = *this;

    for (ZoneType& zn : zones)
    {
        zn.clearAddressing();
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clearPrimitives()
{
    PtrList<ZoneType>& zones = *this;

    for (ZoneType& zn : zones)
    {
        zn.clearPrimitives();
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::clear()
{
    clearAddressing();
    PtrList<ZoneType>::clear();
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::checkDefinition
(
    const bool report
) const
{
    bool hasError = false;

    const PtrList<ZoneType>& zones = *this;

    for (const ZoneType& zn : zones)
    {
        hasError |= zn.checkDefinition(report);
    }

    return hasError;
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::checkParallelSync
(
    const bool report
) const
{
    if (!UPstream::parRun())
    {
        return false;
    }

    const PtrList<ZoneType>& zones = *this;

    bool hasError = false;

    const wordList localNames(this->names());
    const wordList localTypes(this->types());

    // Check and report error(s) on master
    // - don't need indexing on master itself

    const globalIndex procAddr
    (
        globalIndex::gatherNonLocal{},
        localNames.size()
    );

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
                Info<< " ***Inconsistent zones across processors, "
                       "processor 0 has zone names:" << localNames
                    << " zone types:" << localTypes
                    << " processor " << proci
                    << " has zone names:" << procNames
                    << " zone types:" << procTypes
                    << endl;
            }
        }
    }

    Pstream::broadcast(hasError);

    // Check local contents
    if (!hasError)
    {
        for (const ZoneType& zn : zones)
        {
            if (zn.checkParallelSync(false))
            {
                hasError = true;

                if (debug || (report && UPstream::master()))
                {
                    Info<< " ***Zone " << zn.name()
                        << " of type " << zn.type()
                        << " is not correctly synchronised"
                        << " across coupled boundaries."
                        << " (coupled faces are either not both"
                        << " present in set or have same flipmap)" << endl;
                }
            }
        }
    }

    return hasError;
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::movePoints(const pointField& pts)
{
    PtrList<ZoneType>& zones = *this;

    for (ZoneType& zn : zones)
    {
        zn.movePoints(pts);
    }
}


template<class ZoneType, class MeshType>
void Foam::ZoneMesh<ZoneType, MeshType>::updateMetaData()
{
    wordList zoneNames(this->names());
    if (zoneNames.empty())
    {
        this->removeMetaData();
    }
    else
    {
        dictionary& meta = this->getMetaData();
        meta.set("names", zoneNames);
    }
}


template<class ZoneType, class MeshType>
bool Foam::ZoneMesh<ZoneType, MeshType>::writeData(Ostream& os) const
{
    os  << *this;
    return os.good();
}


// * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
const ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator[]
(
    const word& zoneName
) const
{
    const label zonei = findZoneID(zoneName);

    if (zonei < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zonei);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator[]
(
    const word& zoneName
)
{
    const label zonei = findZoneID(zoneName);

    if (zonei < 0)
    {
        FatalErrorInFunction
            << "Zone named " << zoneName << " not found." << nl
            << "Available zone names: " << names() << endl
            << abort(FatalError);
    }

    return operator[](zonei);
}


template<class ZoneType, class MeshType>
ZoneType& Foam::ZoneMesh<ZoneType, MeshType>::operator()
(
    const word& zoneName,
    const bool verbose
)
{
    ZoneType* ptr = findZone(zoneName);

    const bool existing = bool(ptr);

    if (!ptr)
    {
        ptr = new ZoneType(zoneName, this->size(), *this);
        this->push_back(ptr);
    }

    if (verbose)
    {
        Info<< ZoneType::typeName << ' ' << zoneName
            << " (" << (existing ? "existing" : "new")
            << " at index " << ptr->index() << ')'
            << endl;
    }

    return *ptr;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ZoneType, class MeshType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const ZoneMesh<ZoneType, MeshType>& zones
)
{
    const label sz = zones.size();

    if (sz)
    {
        os  << sz << nl << token::BEGIN_LIST;

        for (label i=0; i < sz; ++i)
        {
            zones[i].writeDict(os);
        }

        os  << token::END_LIST;
    }
    else
    {
        os  << sz << token::BEGIN_LIST << token::END_LIST;
    }

    return os;
}


// ************************************************************************* //
