/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2019 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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
#include "predicates.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(objectRegistry, 0);
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Templated implementation for erase() with iterator range.
// Prefer not to expose directly.
template<class InputIter>
static label eraseImpl(objectRegistry& obr, InputIter first, InputIter last)
{
    label changed = 0;

    for
    (
        const label nTotal = obr.size();
        changed < nTotal && first != last; // Terminate early
        ++first
    )
    {
        if (obr.erase(*first))
        {
            ++changed;
        }
    }

    return changed;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::objectRegistry::parentNotTime() const noexcept
{
    return (&parent_ != static_cast<const objectRegistry*>(&time_));
}


// * * * * * * * * * * * * * * * * Constructors *  * * * * * * * * * * * * * //

Foam::objectRegistry::objectRegistry
(
    const Time& t,
    const label initialCapacity
)
:
    regIOobject
    (
        IOobject
        (
            word::validate(t.caseName()),
            t.path(),
            t,
            IOobjectOption::NO_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        true    // to flag that this is the top-level regIOobject
    ),
    HashTable<regIOobject*>(initialCapacity),
    time_(t),
    parent_(t),
    dbDir_(name()),
    event_(1),
    cacheTemporaryObjectsActive_(false),
    cacheTemporaryObjects_(0),
    temporaryObjects_(0)
{}


Foam::objectRegistry::objectRegistry
(
    const IOobject& io,
    const label initialCapacity
)
:
    regIOobject(io),
    HashTable<regIOobject*>(initialCapacity),
    time_(io.time()),
    parent_(io.db()),
    dbDir_(parent_.dbDir()/local()/name()),
    event_(1),
    cacheTemporaryObjectsActive_(false),
    cacheTemporaryObjects_(0),
    temporaryObjects_(0)
{
    writeOpt(IOobjectOption::AUTO_WRITE);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::objectRegistry::~objectRegistry()
{
    objectRegistry::clear();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::objectRegistry::isTimeDb() const noexcept
{
    return (this == &static_cast<const objectRegistry&>(time_));
}


Foam::HashTable<Foam::wordHashSet> Foam::objectRegistry::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::IOobject Foam::objectRegistry::newIOobject
(
    const word& name,
    IOobjectOption ioOpt
) const
{
    return IOobject
    (
        name,
        time().timeName(),  // instance
        *this,
        ioOpt
    );
}


Foam::IOobject Foam::objectRegistry::newIOobject
(
    const word& name,
    IOobjectOption::readOption rOpt,
    IOobjectOption::writeOption wOpt,
    IOobjectOption::registerOption regOpt
) const
{
    return IOobject
    (
        name,
        time().timeName(),  // instance
        *this,
        IOobjectOption(rOpt, wOpt, regOpt)
    );
}


// FUTURE?
//
// Foam::IOobject Foam::objectRegistry::newIOobject_constant
// (
//     const word& name,
//     IOobjectOption::readOption rOpt,
//     IOobjectOption::writeOption wOpt,
//     IOobjectOption::registerOption regOpt
// ) const
// {
//     return IOobject
//     (
//         name,
//         time().constant(),  // instance
//         *this,
//         IOobjectOption(rOpt, wOpt, regOpt)
//     );
// }


Foam::label Foam::objectRegistry::count(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return count(static_cast<word>(clsName));
}


Foam::wordList Foam::objectRegistry::names() const
{
    return HashTable<regIOobject*>::toc();
}


Foam::wordList Foam::objectRegistry::sortedNames() const
{
    return HashTable<regIOobject*>::sortedToc();
}


Foam::wordList Foam::objectRegistry::names(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return names(static_cast<word>(clsName));
}


Foam::wordList Foam::objectRegistry::sortedNames(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return sortedNames(static_cast<word>(clsName));
}


const Foam::objectRegistry& Foam::objectRegistry::subRegistry
(
    const word& name,
    const bool forceCreate,
    const bool recursive
) const
{
    if (forceCreate && !foundObject<objectRegistry>(name, recursive))
    {
        objectRegistry* subObr = new objectRegistry
        (
            IOobject
            (
                name,
                time().constant(),
                *this,
                IOobjectOption::NO_READ,
                IOobjectOption::NO_WRITE,
                IOobjectOption::REGISTER
            )
        );
        subObr->store();
    }

    return lookupObject<objectRegistry>(name, recursive);
}


Foam::label Foam::objectRegistry::getEvent() const
{
    label curEvent = event_++;

    if (event_ == labelMax)
    {
        if (objectRegistry::debug)
        {
            WarningInFunction
                << "Event counter has overflowed. "
                << "Resetting counter on all dependent objects." << nl
                << "This might cause extra evaluations." << endl;
        }

        // Reset event counter
        curEvent = 1;
        event_ = 2;

        // No need to reset dependent objects; overflow is now handled
        // in regIOobject::upToDate
    }

    return curEvent;
}


bool Foam::objectRegistry::checkIn(regIOobject* io) const
{
    if (!io) return false;

    objectRegistry& obr = const_cast<objectRegistry&>(*this);

    if (objectRegistry::debug)
    {
        Pout<< "objectRegistry::checkIn : "
            << name() << " : checking in " << io->name()
            << " of type " << io->type()
            << endl;
    }

    // Delete cached object if it has the same name as io, and in the
    // cacheTemporaryObjects list
    if (!cacheTemporaryObjects_.empty())
    {
        auto cacheIter = cacheTemporaryObjects_.find(io->name());

        if (cacheIter.good())
        {
            iterator iter = obr.find(io->name());

            if
            (
                iter.good()
             && iter.val() != io
             && iter.val()->ownedByRegistry()
            )
            {
                if (objectRegistry::debug)
                {
                    Pout<< "objectRegistry::checkIn : "
                        << name() << " : deleting cached object "
                        << io->name() << endl;
                }

                cacheIter.val().first() = false;
                deleteCachedObject(iter.val());
            }
        }
    }


    bool ok = obr.insert(io->name(), io);

    if (!ok && objectRegistry::debug)
    {
        WarningInFunction
            << name() << " : Attempt to checkIn object with name "
            << io->name() << " which was already checked in"
            << endl;
    }

    return ok;
}


bool Foam::objectRegistry::checkOut(regIOobject* io) const
{
    if (!io) return false;

    objectRegistry& obr = const_cast<objectRegistry&>(*this);

    iterator iter = obr.find(io->name());

    if (iter.good())
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::checkOut : "
                << name() << " : checking out " << io->name()
                << " of type " << io->type()
                << endl;
        }

        if (iter.val() != io)
        {
            if (objectRegistry::debug)
            {
                WarningInFunction
                    << name() << " : Attempt to checkOut copy of "
                    << iter.key()
                    << endl;
            }

            return false;
        }

        return obr.erase(iter);
    }


    if (objectRegistry::debug)
    {
        Pout<< "objectRegistry::checkOut : "
            << name() << " : could not find " << io->name() << " in registry"
            << endl;
    }

    return false;
}


bool Foam::objectRegistry::checkIn(regIOobject& io) const
{
    return checkIn(&io);
}


bool Foam::objectRegistry::checkOut(regIOobject& io) const
{
    return checkOut(&io);
}


bool Foam::objectRegistry::checkOut(const word& key) const
{
    return const_cast<objectRegistry&>(*this).erase(key);
}


void Foam::objectRegistry::clear()
{
    // Free anything owned by the registry, but first unset both
    // 'ownedByRegistry' and 'registered' flags to ensure that the
    // regIOobject destructor will not affect the registry

    for (iterator iter = begin(); iter != end(); ++iter)
    {
        regIOobject* ptr = iter.val();

        if (ptr && ptr->ownedByRegistry())
        {
            if (objectRegistry::debug)
            {
                Pout<< "objectRegistry::clear : " << ptr->name() << nl;
            }

            ptr->release(true);     // Relinquish ownership and registration
            delete ptr;             // Delete also clears fileHandler watches
        }
    }

    HashTable<regIOobject*>::clear();
}


void Foam::objectRegistry::clearStorage()
{
    objectRegistry::clear();
    HashTable<regIOobject*>::clearStorage();
}


bool Foam::objectRegistry::erase(const iterator& iter)
{
    // Remove from registry - see notes in objectRegistry::clear()

    if (iter.good())
    {
        regIOobject* ptr = const_cast<iterator&>(iter).val();

        const bool ok = HashTable<regIOobject*>::erase(iter);

        if (ptr && ptr->ownedByRegistry())
        {
            ptr->release(true);     // Relinquish ownership and registration
            delete ptr;             // Delete also clears fileHandler watches
        }

        return ok;
    }

    return false;
}


bool Foam::objectRegistry::erase(const word& key)
{
    return erase(find(key));
}


Foam::label Foam::objectRegistry::erase(std::initializer_list<word> keys)
{
    return eraseImpl(*this, keys.begin(), keys.end());
}


Foam::label Foam::objectRegistry::erase(const UList<word>& keys)
{
    return eraseImpl(*this, keys.begin(), keys.end());
}


void Foam::objectRegistry::rename(const word& newName)
{
    regIOobject::rename(newName);

    // Adjust dbDir_ as well
    const auto i = dbDir_.rfind('/');

    if (i == string::npos)
    {
        dbDir_ = newName;
    }
    else
    {
        dbDir_.replace(i+1, string::npos, newName);
    }
}


const Foam::regIOobject* Foam::objectRegistry::cfindIOobject
(
    const word& name,
    const bool recursive
) const
{
    const_iterator iter = cfind(name);

    if (iter.good())
    {
        return iter.val();
    }
    else if (recursive && !this->isTimeDb())
    {
        return parent_.cfindIOobject(name, recursive);
    }

    return nullptr;
}


bool Foam::objectRegistry::contains
(
    const word& name,
    const bool recursive
) const
{
    return cfindIOobject(name, recursive);
}


bool Foam::objectRegistry::modified() const
{
    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if (iter.val()->modified())
        {
            return true;
        }
    }

    return false;
}


void Foam::objectRegistry::readModifiedObjects()
{
    for (iterator iter = begin(); iter != end(); ++iter)
    {
        if (objectRegistry::debug)
        {
            Pout<< "objectRegistry::readModifiedObjects() : "
                << name() << " : Considering reading object "
                << iter.key() << endl;
        }

        iter.val()->readIfModified();
    }
}


bool Foam::objectRegistry::readIfModified()
{
    readModifiedObjects();
    return true;
}


bool Foam::objectRegistry::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = true;

    for (const_iterator iter = cbegin(); iter != cend(); ++iter)
    {
        if (objectRegistry::debug)
        {
            const regIOobject& obj = *iter.val();

            Pout<< "objectRegistry::write() : "
                << name() << " : Considering writing object "
                << iter.key() << " of type "
                << obj.type() << " with writeOpt "
                << static_cast<int>(obj.writeOpt())
                << " to file " << obj.objectRelPath() << endl;
        }

        if (iter.val()->writeOpt() != IOobjectOption::NO_WRITE)
        {
            ok = iter.val()->writeObject(streamOpt, writeOnProc) && ok;
        }
    }

    return ok;
}


// ************************************************************************* //
