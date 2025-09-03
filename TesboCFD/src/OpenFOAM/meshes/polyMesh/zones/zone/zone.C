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

#include "zone.H"
#include "dictionary.H"
#include "HashSet.H"
#include "IOstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(zone);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zone::zone()
:
    zoneIdentifier(),
    labelList(),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone(const word& name, const label index)
:
    zoneIdentifier(name, index),
    labelList(),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const word& name,
    const labelUList& addr,
    const label index
)
:
    zone(name, index)
{
    labelList::operator=(addr);
}


Foam::zone::zone
(
    const word& name,
    labelList&& addr,
    const label index
)
:
    zone(name, index)
{
    labelList::transfer(addr);
}


Foam::zone::zone
(
    const word& name,
    const dictionary& dict,
    const word& labelsName,
    const label index
)
:
    zoneIdentifier(name, dict, index),
    labelList(dict.get<labelList>(labelsName)),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& originalZone,
    const labelUList& addr,
    const label newIndex
)
:
    zoneIdentifier(originalZone, newIndex),
    labelList(addr),
    lookupMapPtr_(nullptr)
{}


Foam::zone::zone
(
    const zone& originalZone,
    labelList&& addr,
    const label newIndex
)
:
    zoneIdentifier(originalZone, newIndex),
    labelList(std::move(addr)),
    lookupMapPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Map<Foam::label>& Foam::zone::lookupMap() const
{
    if (!lookupMapPtr_)
    {
        const labelList& addr = *this;

        lookupMapPtr_.reset(new Map<label>(2*addr.size()));
        auto& map = *lookupMapPtr_;

        for (const label id : addr)
        {
            map.insert(id, map.size());
        }
    }

    return *lookupMapPtr_;
}


Foam::label Foam::zone::localID(const label globalID) const
{
    return lookupMap().lookup(globalID, -1);
}


void Foam::zone::clearAddressing()
{
    lookupMapPtr_.reset(nullptr);
}


void Foam::zone::clearPrimitives()
{
    static_cast<labelList&>(*this).clear();
}


bool Foam::zone::checkDefinition(const label maxSize, const bool report) const
{
    const labelList& addr = *this;

    bool hasError = false;

    // To check for duplicate entries
    labelHashSet elems(2*size());

    for (const label id : addr)
    {
        if (id < 0 || id >= maxSize)
        {
            hasError = true;

            if (report)
            {
                SeriousErrorInFunction
                    << "Zone " << this->name()
                    << " contains invalid index label " << id << nl
                    << "Valid index labels are 0.."
                    << maxSize-1 << endl;
            }
            else
            {
                // w/o report - can stop checking now
                break;
            }
        }
        else if (!elems.insert(id))
        {
            if (report)
            {
                WarningInFunction
                    << "Zone " << this->name()
                    << " contains duplicate index label " << id << endl;
            }
        }
    }

    return hasError;
}


void Foam::zone::write(Ostream& os) const
{
    os  << nl << this->name()
        << nl << static_cast<const labelList&>(*this);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const zone& zn)
{
    zn.write(os);
    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
