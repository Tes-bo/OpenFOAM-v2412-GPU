/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "boundaryRegion.H"
#include "IOMap.H"
#include "OFstream.H"
#include "predicates.H"

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
            lookupName = "boundaryRegion_" + Foam::name(iter.key());
        }

        if (nameMatcher(lookupName))
        {
            output.emplace(iter.key(), std::move(lookupName));
        }
    }

    return output;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundaryRegion::boundaryRegion
(
    const objectRegistry& obr,
    const word& name,
    const fileName& instance
)
{
    readDict(obr, name, instance);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::label Foam::boundaryRegion::maxIndex() const
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


Foam::label Foam::boundaryRegion::push_back(const dictionary& dict)
{
    label maxId = this->maxIndex();

    insert(++maxId, dict);
    return maxId;
}


Foam::Map<Foam::word> Foam::boundaryRegion::names() const
{
    return names_impl(*this, predicates::always{});
}


Foam::Map<Foam::word> Foam::boundaryRegion::names(const wordRes& patterns) const
{
    return names_impl(*this, patterns);
}


Foam::Map<Foam::word> Foam::boundaryRegion::boundaryTypes() const
{
    Map<word> output;
    output.reserve(size());

    forAllConstIters(*this, iter)
    {
        word lookupType;
        if (!iter().readIfPresent("BoundaryType", lookupType))
        {
            lookupType = "patch";
        }

        output.emplace(iter.key(), std::move(lookupType));
    }

    return output;
}


Foam::word Foam::boundaryRegion::name(const label id) const
{
    word lookupName;

    const auto iter = cfind(id);
    if (iter.good())
    {
        iter.val().readIfPresent("Label", lookupName);
    }

    if (lookupName.empty() && id >= 0)
    {
        lookupName = "boundaryRegion_" + Foam::name(id);
    }

    return lookupName;
}


Foam::label Foam::boundaryRegion::findIndex(const word& name) const
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


Foam::word Foam::boundaryRegion::boundaryType(const word& name) const
{
    word lookupType("patch");

    const label id = this->findIndex(name);
    if (id >= 0)
    {
        operator[](id).readIfPresent("BoundaryType", lookupType);
    }

    return lookupType;
}


void Foam::boundaryRegion::readDict
(
    const objectRegistry& obr,
    const word& name,
    const fileName& instance
)
{
    Map<dictionary>::clear();

    // Read constant/dictName
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
    }
    else
    {
        Info<< "no constant/boundaryRegion information available" << endl;
    }
}


void Foam::boundaryRegion::writeDict
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

void Foam::boundaryRegion::operator=(const boundaryRegion& rhs)
{
    Map<dictionary>::operator=(rhs);
}


void Foam::boundaryRegion::operator=(const Map<dictionary>& rhs)
{
    Map<dictionary>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

void Foam::boundaryRegion::rename(const dictionary& mapDict)
{
    if (mapDict.empty())
    {
        return;
    }

    // Use 1st pass to collect all the regions to be changed
    // and 2nd pass to relabel regions.
    // This avoid re-matching any renamed regions

    Map<word> mapping;
    mapping.reserve(mapDict.size());

    for (const entry& dEntry : mapDict)
    {
        const word oldName(dEntry.stream());

        const label id = this->findIndex(oldName);
        if (id >= 0)
        {
            mapping.insert(id, dEntry.keyword());
        }
    }

    forAllConstIters(mapping, iter)
    {
        const word& newName = iter.val();

        dictionary& dict = operator[](iter.key());

        word oldName;
        if (!dict.readIfPresent("Label", oldName))
        {
            oldName = "boundaryRegion_" + Foam::name(iter.key());
        }

        Info<< "rename patch: " << newName << " <- " << oldName << nl;

        dict.set("Label", newName);
    }
}


// ************************************************************************* //
