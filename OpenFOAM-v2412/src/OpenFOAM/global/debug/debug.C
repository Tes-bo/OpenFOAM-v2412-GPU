/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Description
    Class for handling debugging switches.

Note
    Included by global/globals.C

\*---------------------------------------------------------------------------*/

#include "debug.H"
#include "dictionary.H"
#include "IFstream.H"
#include "etcFiles.H"
#include "Ostream.H"
#include "demandDrivenData.H"
#include "simpleObjectRegistry.H"
#include "IOobject.H"
#include "HashSet.H"
#include "nullObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace debug
{

//! \cond ignoreDocumentation
//- Skip documentation : local scope only

dictionary* controlDictPtr_(nullptr);
dictionary* debugSwitchesPtr_(nullptr);
dictionary* infoSwitchesPtr_(nullptr);
dictionary* optimisationSwitchesPtr_(nullptr);

// Debug switch read and write callback tables.
simpleObjectRegistry* debugObjectsPtr_(nullptr);
simpleObjectRegistry* infoObjectsPtr_(nullptr);
simpleObjectRegistry* optimisationObjectsPtr_(nullptr);
simpleObjectRegistry* dimensionSetObjectsPtr_(nullptr);
simpleObjectRegistry* dimensionedConstantObjectsPtr_(nullptr);


// To ensure controlDictPtr_ is deleted at the end of the run
struct deleteControlDictPtr
{
    ~deleteControlDictPtr()
    {
        deleteDemandDrivenData(debugObjectsPtr_);
        deleteDemandDrivenData(infoObjectsPtr_);
        deleteDemandDrivenData(optimisationObjectsPtr_);
        deleteDemandDrivenData(dimensionSetObjectsPtr_);
        deleteDemandDrivenData(dimensionedConstantObjectsPtr_);

        debugSwitchesPtr_ = nullptr;
        infoSwitchesPtr_ = nullptr;
        optimisationSwitchesPtr_ = nullptr;
        deleteDemandDrivenData(controlDictPtr_);
    }
};

deleteControlDictPtr deleteControlDictPtr_;
//! \endcond


} // End namespace debug
} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Like dictionary getOrAdd with LITERAL, but circumventing
// writeOptionalEntries to avoid extremely noisy output
template<class T>
static inline T getOrAdd
(
    dictionary& dict,
    const char* name,
    const T deflt
)
{
    const entry* eptr = dict.findEntry(name, keyType::LITERAL);

    if (eptr)
    {
        return eptr->get<T>();
    }

    dict.add(new primitiveEntry(name, deflt));
    return deflt;
}


// Append object to a registry
static inline void appendNamedEntry
(
    simpleObjectRegistry& obr,
    const char* name,
    simpleRegIOobject* obj
)
{
    simpleObjectRegistryEntry* ptr = obr.find(name);
    if (ptr)
    {
        ptr->append(obj);
    }
    else
    {
        obr.append(name, new simpleObjectRegistryEntry(obj));
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::dictionary& Foam::debug::controlDict()
{
    if (!controlDictPtr_)
    {
        string controlDictString(Foam::getEnv("FOAM_CONTROLDICT"));
        if (!controlDictString.empty())
        {
            // Read from environment
            IStringStream is(controlDictString);
            controlDictPtr_ = new dictionary(is);
        }
        else
        {
            fileNameList controlDictFiles = findEtcFiles("controlDict", true);
            controlDictPtr_ = new dictionary();
            forAllReverse(controlDictFiles, i)
            {
                IFstream is(controlDictFiles[i]);

                if (!is.good())
                {
                    SafeFatalIOErrorInFunction
                    (
                        is,
                        "Cannot open controlDict"
                    );
                }
                controlDictPtr_->merge(dictionary(is));
            }
        }
    }

    return *controlDictPtr_;
}


Foam::dictionary& Foam::debug::switchSet
(
    const char* subDictName,
    dictionary*& subDictPtr
)
{
    if (!subDictPtr)
    {
        subDictPtr = controlDict().findDict(subDictName, keyType::LITERAL);

        if (!subDictPtr)
        {
            std::cerr
                << "debug::switchSet(const char*, dictionary*&):\n"
                << "    Cannot find " <<  subDictName << " in dictionary "
                << controlDict().name().c_str()
                << std::endl << std::endl;

            std::exit(1);
        }
    }

    return *subDictPtr;
}


Foam::dictionary& Foam::debug::debugSwitches()
{
    return switchSet("DebugSwitches", debugSwitchesPtr_);
}


Foam::dictionary& Foam::debug::infoSwitches()
{
    return switchSet("InfoSwitches", infoSwitchesPtr_);
}


Foam::dictionary& Foam::debug::optimisationSwitches()
{
    return switchSet("OptimisationSwitches", optimisationSwitchesPtr_);
}


int Foam::debug::debugSwitch(const char* name, const int deflt)
{
    return getOrAdd(debugSwitches(), name, deflt);
}


int Foam::debug::infoSwitch(const char* name, const int deflt)
{
    return getOrAdd(infoSwitches(), name, deflt);
}


int Foam::debug::optimisationSwitch(const char* name, const int deflt)
{
    return getOrAdd(optimisationSwitches(), name, deflt);
}


float Foam::debug::floatOptimisationSwitch(const char* name, const float deflt)
{
    return getOrAdd(optimisationSwitches(), name, deflt);
}


void Foam::debug::addDebugObject(const char* name, simpleRegIOobject* obj)
{
    appendNamedEntry(debugObjects(), name, obj);
}


void Foam::debug::addInfoObject(const char* name, simpleRegIOobject* obj)
{
    appendNamedEntry(infoObjects(), name, obj);
}


void Foam::debug::addOptimisationObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    appendNamedEntry(optimisationObjects(), name, obj);
}


void Foam::debug::addDimensionSetObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    appendNamedEntry(dimensionSetObjects(), name, obj);
}


void Foam::debug::addDimensionedConstantObject
(
    const char* name,
    simpleRegIOobject* obj
)
{
    appendNamedEntry(dimensionedConstantObjects(), name, obj);
}


Foam::simpleObjectRegistry& Foam::debug::debugObjects()
{
    if (!debugObjectsPtr_)
    {
        debugObjectsPtr_ = new simpleObjectRegistry(128);
    }

    return *debugObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::infoObjects()
{
    if (!infoObjectsPtr_)
    {
        infoObjectsPtr_ = new simpleObjectRegistry(128);
    }

    return *infoObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::optimisationObjects()
{
    if (!optimisationObjectsPtr_)
    {
        optimisationObjectsPtr_ = new simpleObjectRegistry(128);
    }

    return *optimisationObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionSetObjects()
{
    if (!dimensionSetObjectsPtr_)
    {
        dimensionSetObjectsPtr_ = new simpleObjectRegistry(128);
    }

    return *dimensionSetObjectsPtr_;
}


Foam::simpleObjectRegistry& Foam::debug::dimensionedConstantObjects()
{
    if (!dimensionedConstantObjectsPtr_)
    {
        dimensionedConstantObjectsPtr_ = new simpleObjectRegistry(128);
    }

    return *dimensionedConstantObjectsPtr_;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Print the switch status
static inline void printStatus
(
    const char * const message,
    const wordList& list
)
{
    // Use writeList with length = -1 to ensure we always have newlines,
    // even for short lists

    Info<< message << nl;
    list.writeList(Info, -1) << nl;
}


// Write the switch names.
//
// Use writeList with -1 for the length to ensure we always have newlines,
// even if the lists are short

static void listSwitches
(
    const wordList& debugSwitches,
    const wordList& infoSwitches,
    const wordList& optSwitches,
    const bool unset
)
{
    IOobject::writeDivider(Info);

    if (unset)
    {
        fileNameList controlDictFiles = findEtcFiles("controlDict", true);
        dictionary controlDict;
        forAllReverse(controlDictFiles, i)
        {
            IFstream is(controlDictFiles[i]);

            controlDict.merge(dictionary(is));
        }

        // HashSet to track switches that have not been set
        wordHashSet hashed;

        // DebugSwitches
        if (notNull(debugSwitches))
        {
            hashed = debugSwitches;
            hashed.unset(controlDict.subDict("DebugSwitches").toc());
            printStatus("Unset DebugSwitches", hashed.sortedToc());
        }

        // InfoSwitches
        if (notNull(infoSwitches))
        {
            hashed = infoSwitches;
            hashed.unset(controlDict.subDict("InfoSwitches").toc());
            printStatus("Unset InfoSwitches", hashed.sortedToc());
        }

        // OptimisationSwitches
        if (notNull(optSwitches))
        {
            hashed = optSwitches;
            hashed.unset(controlDict.subDict("OptimisationSwitches").toc());
            printStatus("Unset OptimisationSwitches", hashed.sortedToc());
        }
    }
    else
    {
        // DebugSwitches
        if (notNull(debugSwitches))
        {
            printStatus("DebugSwitches", debugSwitches);
        }

        // InfoSwitches
        if (notNull(infoSwitches))
        {
            printStatus("InfoSwitches", infoSwitches);
        }

        // OptimisationSwitches
        if (notNull(optSwitches))
        {
            printStatus("OptimisationSwitches", optSwitches);
        }
    }
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::debug::listSwitches(const bool unset)
{
    listSwitches
    (
        debug::debugSwitches().sortedToc(),
        debug::infoSwitches().sortedToc(),
        debug::optimisationSwitches().sortedToc(),
        unset
    );
}


void Foam::debug::listDebugSwitches(const bool unset)
{
    listSwitches
    (
        debug::debugSwitches().sortedToc(),
        wordList::null(),
        wordList::null(),
        unset
    );
}


void Foam::debug::listInfoSwitches(const bool unset)
{
    listSwitches
    (
        wordList::null(),
        debug::infoObjects().sortedToc(),
        wordList::null(),
        unset
    );
}


void Foam::debug::listOptimisationSwitches(const bool unset)
{
    listSwitches
    (
        wordList::null(),
        wordList::null(),
        debug::optimisationSwitches().sortedToc(),
        unset
    );
}


void Foam::debug::listRegisteredSwitches(const bool unset)
{
    listSwitches
    (
        debug::debugObjects().sortedToc(),
        debug::infoObjects().sortedToc(),
        debug::optimisationObjects().sortedToc(),
        unset
    );
}


void Foam::debug::listRegisteredDebugSwitches(const bool unset)
{
    listSwitches
    (
        debug::debugObjects().sortedToc(),
        wordList::null(),
        wordList::null(),
        unset
    );
}


void Foam::debug::listRegisteredInfoSwitches(const bool unset)
{
    listSwitches
    (
        wordList::null(),
        debug::infoObjects().sortedToc(),
        wordList::null(),
        unset
    );
}


void Foam::debug::listRegisteredOptimisationSwitches(const bool unset)
{
    listSwitches
    (
        wordList::null(),
        wordList::null(),
        debug::optimisationObjects().sortedToc(),
        unset
    );
}


// ************************************************************************* //
