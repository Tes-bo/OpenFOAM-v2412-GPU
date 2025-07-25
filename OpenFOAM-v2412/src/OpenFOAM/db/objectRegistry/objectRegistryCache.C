/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::objectRegistry::readCacheTemporaryObjects() const
{
    if (cacheTemporaryObjectsActive_) return;

    const auto* eptr = time_.controlDict().findEntry
    (
        "cacheTemporaryObjects",
        keyType::LITERAL
    );

    if (eptr)
    {
        cacheTemporaryObjectsActive_ = true;

        // Clear old cache?
        // cacheTemporaryObjects_.clear();

        wordList objectNames;

        if (eptr->isDict())
        {
            // Per region (sub-dictionary syntax)
            eptr->dict().readIfPresent(name(), objectNames);
        }
        else
        {
            // All regions
            eptr->readEntry(objectNames);
        }

        for (const word& objName : objectNames)
        {
            cacheTemporaryObjects_.emplace(objName, false, false);
        }

        if (objectRegistry::debug)
        {
            Info<< "objectRegistry::cacheTemporaryObjects : "
                << flatOutput(objectNames) << endl;
        }
    }
}


void Foam::objectRegistry::deleteCachedObject(regIOobject* io) const
{
    if (io)
    {
        io->release();     // Relinquish any ownership by registry
        io->checkOut();
        // Additional safety - not certain this is actually needed...
        io->rename(io->name() + "-Cache");
        delete io;
    }
}


// FUTURE: (currently not needed)
// void Foam::objectRegistry::addTemporaryObject
// (
//     const word& name
// ) const
// {
//     cacheTemporaryObjects_.emplace(name, false, false);
// }


bool Foam::objectRegistry::is_cacheTemporaryObject
(
    const word& name
) const
{
    return cacheTemporaryObjects_.contains(name);
}


bool Foam::objectRegistry::is_cacheTemporaryObject(const regIOobject* io) const
{
    return io && cacheTemporaryObjects_.contains(io->name());
}


bool Foam::objectRegistry::is_cacheTemporaryObject(const regIOobject& io) const
{
    return cacheTemporaryObjects_.contains(io.name());
}


void Foam::objectRegistry::resetCacheTemporaryObject
(
    const regIOobject* io
) const
{
    if (io && !cacheTemporaryObjects_.empty())
    {
        auto iter = cacheTemporaryObjects_.find(io->name());

        // Reset the cached flag
        if (iter.good())
        {
            iter.val().first() = false;
        }
    }
}


void Foam::objectRegistry::resetCacheTemporaryObject
(
    const regIOobject& io
) const
{
    resetCacheTemporaryObject(&io);
}


bool Foam::objectRegistry::checkCacheTemporaryObjects() const
{
    bool enabled = cacheTemporaryObjects_.size();

    forAllConstIters(*this, iter)
    {
        const auto* subObr = dynamic_cast<const objectRegistry*>(iter.val());

        // Protect against re-searching the top-level registry
        if (subObr && subObr != this)
        {
            enabled = subObr->checkCacheTemporaryObjects() || enabled;
        }
    }

    if (enabled)
    {
        OSstream* emitWarning = nullptr;

        forAllIters(cacheTemporaryObjects_, iter)
        {
            if (!iter.val().second())
            {
                if (!emitWarning)
                {
                    emitWarning = &(Foam::Warning.stream());

                    *emitWarning
                        << "objectRegistry '"
                        << name() << "' has missing temporary objects:" << nl;
                }

                *emitWarning<< "    " << iter.key() << nl;
            }
            else
            {
                iter.val().second() = false;
            }
        }

        if (emitWarning)
        {
            *emitWarning
                << "Available temporary objects: "
                << temporaryObjects_.sortedToc() << endl;
        }

        temporaryObjects_.clear();
    }

    return enabled;
}


// ************************************************************************* //
