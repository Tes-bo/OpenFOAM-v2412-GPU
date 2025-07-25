/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

#include "cellTable.H"
#include "IOMap.H"
#include "polyMesh.H"
#include "OFstream.H"
#include "predicates.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

static const char* const defaultMaterial_ = "fluid";


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

template<class MatchPredicate>
static Map<word> names_impl
(
    const Map<dictionary>& input,
    const MatchPredicate& nameMatcher
)
{
    Map<word> output;
    output.reserve(input.size());

    forAllConstIters(input, iter)
    {
        word lookupName;
        if (!iter().readIfPresent("Label", lookupName))
        {
            lookupName = "cellTable_" + Foam::name(iter.key());
        }

        if (nameMatcher(lookupName))
        {
            output.emplace(iter.key(), std::move(lookupName));
        }
    }

    return output;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellTable::addDefaults()
{
    forAllIters(*this, iter)
    {
        if (!iter().found("MaterialType"))
        {
            iter().add("MaterialType", word(defaultMaterial_));
        }
    }
}


void Foam::cellTable::setEntry
(
    const label id,
    const word& key,
    const word& value
)
{
    dictionary dict;
    dict.add(key, value);

    iterator iter = find(id);
    if (iter.good())
    {
        iter().merge(dict);
    }
    else
    {
        insert(id, dict);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cellTable::cellTable
(
    const objectRegistry& obr,
    const word& name,
    const fileName& instance
)
{
    readDict(obr, name, instance);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::cellTable::maxIndex() const
{
    label maxId = -1;
    forAllConstIters(*this, iter)
    {
        if (maxId < iter.key())
        {
            maxId = iter.key();
        }
    }

    return maxId;
}


Foam::label Foam::cellTable::push_back(const dictionary& dict)
{
    label maxId = this->maxIndex();

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::cellTable::names() const
{
    return names_impl(*this, predicates::always{});
}


Foam::Map<Foam::word> Foam::cellTable::names(const wordRes& patterns) const
{
    return names_impl(*this, patterns);
}


Foam::word Foam::cellTable::name(const label id) const
{
    word lookupName;

    const auto iter = cfind(id);
    if (iter.good())
    {
        iter.val().readIfPresent("Label", lookupName);
    }

    if (lookupName.empty() && id >= 0)
    {
        lookupName = "cellTable_" + Foam::name(id);
    }

    return lookupName;
}


Foam::label Foam::cellTable::findIndex(const word& name) const
{
    if (name.empty())
    {
        return -1;
    }

    forAllConstIters(*this, iter)
    {
        const auto& dict = iter.val();

        word lookupName;
        if (dict.readIfPresent("Label", lookupName) && (lookupName == name))
        {
            return iter.key();
        }
    }

    return -1;
}


Foam::Map<Foam::word> Foam::cellTable::materialTypes() const
{
    Map<word> output;
    output.reserve(size());

    forAllConstIters(*this, iter)
    {
        word lookupType;
        if (!iter().readIfPresent("MaterialType", lookupType))
        {
            lookupType = defaultMaterial_;
        }

        output.emplace(iter.key(), std::move(lookupType));
    }

    return output;
}


Foam::Map<Foam::word> Foam::cellTable::selectType(const word& matl) const
{
    Map<word> output;
    output.reserve(size());

    forAllConstIters(*this, iter)
    {
        const dictionary& dict = iter.val();

        if
        (
            matl
         == dict.getOrDefault<word>("MaterialType", defaultMaterial_)
        )
        {
            word lookupName;
            if (dict.readIfPresent("Label", lookupName))
            {
                lookupName = "cellTable_" + Foam::name(iter.key());
            }

            output.emplace(iter.key(), std::move(lookupName));
        }
    }

    return output;
}


Foam::Map<Foam::word> Foam::cellTable::fluids() const
{
    return selectType("fluid");
}


Foam::Map<Foam::word> Foam::cellTable::solids() const
{
    return selectType("solid");
}


Foam::Map<Foam::word> Foam::cellTable::shells() const
{
    return selectType("shell");
}


void Foam::cellTable::setMaterial(const label id, const word& matlType)
{
    setEntry(id, "MaterialType", matlType);
}


void Foam::cellTable::setName(const label id, const word& name)
{
    setEntry(id, "Label", name);
}


void Foam::cellTable::setName(const label id)
{
    iterator iter = find(id);

    if (!iter.good() || !iter().found("Label"))
    {
        setName(id, "cellTable_" + Foam::name(id));
    }
}


void Foam::cellTable::readDict
(
    const objectRegistry& obr,
    const word& name,
    const fileName& instance
)
{
    clear();

    // read constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            obr,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    if (ioObj.headerOk())
    {
        *this = ioObj;
        addDefaults();
    }
    else
    {
        Info<< "no constant/cellTable information available" << endl;
    }
}


void Foam::cellTable::writeDict
(
    const objectRegistry& obr,
    const word& name,
    const fileName& instance
) const
{
    // write constant/dictName
    IOMap<dictionary> ioObj
    (
        IOobject
        (
            name,
            instance,
            obr,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    ioObj.note() =
        "persistent data for third-party mesh <-> OpenFOAM translation";

    Info<< "Writing " << ioObj.name() << " to "
        << ioObj.objectRelPath() << endl;

    OFstream os(ioObj.objectPath());
    ioObj.writeHeader(os);
    os << *this;

    IOobject::writeEndDivider(os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::cellTable::operator=(const cellTable& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
    addDefaults();
}


void Foam::cellTable::operator=(const polyMesh& mesh)
{
    Map<dictionary> zoneDict;

    // create cellTableId and cellTable based on cellZones
    label nZoneCells = 0;

    wordList zoneNames = mesh.cellZones().names();
    label unZonedType = zoneNames.size() + 1;

    // do cell zones
    forAll(mesh.cellZones(), zoneI)
    {
        const cellZone& cZone = mesh.cellZones()[zoneI];
        nZoneCells += cZone.size();

        dictionary dict;
        dict.add("Label", zoneNames[zoneI]);
        zoneDict.insert(zoneI + 1, dict);
    }

    // collect unzoned cells
    // special case: no zones at all - do entire mesh
    if (nZoneCells == 0)
    {
        zoneDict.clear();
        unZonedType = 1;
    }

    if (mesh.nCells() > nZoneCells)
    {
        zoneDict.insert
        (
            unZonedType,
            dictionary(IStringStream("Label cells;")())
        );
    }

    Map<dictionary>::operator=(zoneDict);
    addDefaults();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::cellTable::addCellZones
(
    polyMesh& mesh,
    const labelList& tableIds
) const
{
    // From cellTable ID => zone index
    Map<label> typeToZone;

    // Name per zone (not cellTableID)
    wordList zoneNames;

    {
        // The cellTable ID => zone name
        Map<word> namesLookup = this->names();

        zoneNames.resize(namesLookup.size());
        typeToZone.reserve(namesLookup.size());

        // Linear indexing
        label zonei = 0;
        for (const label id : namesLookup.sortedToc())
        {
            typeToZone(id) = zonei;
            zoneNames[zonei] = namesLookup[id];
            ++zonei;
        }
    }


    List<DynamicList<label>> zoneCells(zoneNames.size());

    forAll(tableIds, celli)
    {
        label zonei = typeToZone.lookup(tableIds[celli], -1);
        if (zonei >= 0)
        {
            zoneCells[zonei].push_back(celli);
        }
    }

    // Track which zones were actually used
    DynamicList<label> zoneUsed(zoneCells.size());

    forAll(zoneCells, zonei)
    {
        zoneCells[zonei].shrink();
        if (!zoneCells[zonei].empty())
        {
            zoneUsed.push_back(zonei);
        }
    }

    const label nZonesUsed = zoneUsed.size();

    cellZoneMesh& czMesh = mesh.cellZones();

    czMesh.clear();
    if (nZonesUsed <= 1)
    {
        Info<< "cellZones not used" << endl;
        return;
    }
    czMesh.resize(nZonesUsed);

    forAll(zoneUsed, zonei)
    {
        const label origZonei = zoneUsed[zonei];

        Info<< "cellZone " << zonei
            << " (size: "  << zoneCells[origZonei].size()
            << ") name: "  << zoneNames[origZonei] << endl;

        czMesh.set
        (
            zonei,
            new cellZone
            (
                zoneNames[origZonei],
                zoneCells[origZonei],
                zonei,
                czMesh
            )
        );
    }
    czMesh.writeOpt(IOobject::AUTO_WRITE);
}


void Foam::cellTable::combine(const dictionary& mapDict, labelList& tableIds)
{
    if (mapDict.empty())
    {
        return;
    }

    Map<word> origNames(this->names());
    labelList mapping(identity(this->maxIndex() + 1));

    bool remap = false;
    for (const entry& dEntry : mapDict)
    {
        wordRes patterns(dEntry.stream());

        // find all matches
        Map<word> matches;
        forAllConstIters(origNames, namesIter)
        {
            if (patterns.match(namesIter()))
            {
                matches.insert(namesIter.key(), namesIter());
            }
        }

        if (matches.size())
        {
            label targetId = this->findIndex(dEntry.keyword());

            Info<< "combine cellTable: " << dEntry.keyword();
            if (targetId < 0)
            {
                // not found - reuse 1st element but with different name
                targetId = min(matches.toc());
                operator[](targetId).set("Label", dEntry.keyword());

                Info<< " = (";
            }
            else
            {
                Info<< " += (";
            }


            // the mapping and name for targetId is already okay
            matches.erase(targetId);
            origNames.erase(targetId);

            // remove matched names, leaving targetId on 'this'
            this->erase(matches);
            origNames.erase(matches);

            forAllConstIters(matches, matchIter)
            {
                mapping[matchIter.key()] = targetId;
                Info<< " " << matchIter();
            }
            Info<< " )" << endl;

            remap = true;
        }
    }

    if (remap)
    {
        inplaceRenumber(mapping, tableIds);
    }
}

// ************************************************************************* //
