/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "IOobjectList.H"
#include "Time.H"
#include "predicates.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::IOobjectList::checkObjectOrder
(
    const UPtrList<const IOobject>& objs,
    bool syncPar
)
{
    if (syncPar && UPstream::is_parallel())
    {
        wordList objectNames(objs.size());

        auto iter = objectNames.begin();

        for (const IOobject& io : objs)
        {
            *iter = io.name();   // nameOp<IOobject>()
            ++iter;
        }

        checkNameOrder(objectNames, syncPar);
    }
}


void Foam::IOobjectList::checkNameOrder
(
    const wordList& objectNames,
    bool syncPar
)
{
    if (syncPar && UPstream::is_parallel())
    {
        wordList masterNames;
        if (UPstream::master())
        {
            masterNames = objectNames;
        }
        Pstream::broadcast(masterNames);

        if (!UPstream::master() && (objectNames != masterNames))
        {
            FatalErrorInFunction
                << "Objects not synchronised across processors." << nl
                << "Master has " << flatOutput(masterNames) << nl
                << "Processor " << UPstream::myProcNo()
                << " has " << flatOutput(objectNames) << endl
                << exit(FatalError);
        }
    }
}


void Foam::IOobjectList::syncNames(wordList& objNames)
{
    // Synchronize names
    Pstream::combineReduce(objNames, ListOps::uniqueEqOp<word>());
    Foam::sort(objNames);  // Consistent order
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobjectList::IOobjectList
(
    const objectRegistry& db,
    const fileName& instance,
    const fileName& local,
    IOobjectOption ioOpt
)
:
    HashPtrTable<IOobject>()
{
    word newInstance;
    fileNameList objNames = fileHandler().readObjects
    (
        db,
        instance,
        local,
        newInstance
    );

    for (const auto& objName : objNames)
    {
        auto objectPtr = autoPtr<IOobject>::New
        (
            objName,
            newInstance,
            local,
            db,
            ioOpt
        );

        bool ok = false;
        const bool oldThrowingIOerr = FatalIOError.throwing(true);

        try
        {
            // Use object with local scope and current instance (no searching)
            ok = objectPtr->typeHeaderOk<regIOobject>(false, false);
        }
        catch (const Foam::IOerror& err)
        {
            Warning << err << nl << endl;
        }

        FatalIOError.throwing(oldThrowingIOerr);

        if (ok)
        {
            insert(objectPtr->name(), objectPtr);
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::IOobject* Foam::IOobjectList::cfindObject
(
    const word& objName
) const
{
    // Like HashPtrTable::get(), or lookup() with a nullptr
    const IOobject* io = nullptr;

    const const_iterator iter(cfind(objName));
    if (iter.good())
    {
        io = iter.val();
    }

    if (IOobject::debug)
    {
        if (io)
        {
            InfoInFunction << "Found " << objName << endl;
        }
        else
        {
            InfoInFunction << "Could not find " << objName << endl;
        }
    }

    return io;
}


const Foam::IOobject* Foam::IOobjectList::findObject
(
    const word& objName
) const
{
    return cfindObject(objName);
}


Foam::IOobject* Foam::IOobjectList::findObject(const word& objName)
{
    return const_cast<IOobject*>(cfindObject(objName));
}


Foam::IOobject* Foam::IOobjectList::getObject(const word& objName) const
{
    return const_cast<IOobject*>(cfindObject(objName));
}


Foam::IOobjectList Foam::IOobjectList::lookupClass(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return lookupClass(static_cast<word>(clsName));
}


Foam::HashTable<Foam::wordHashSet> Foam::IOobjectList::classes() const
{
    return classesImpl(*this, predicates::always());
}


Foam::label Foam::IOobjectList::count(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return count(static_cast<word>(clsName));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::wordList Foam::IOobjectList::names() const
{
    return HashPtrTable<IOobject>::toc();
}


Foam::wordList Foam::IOobjectList::names(const bool syncPar) const
{
    return sortedNames(syncPar);
}


Foam::wordList Foam::IOobjectList::names(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return names(static_cast<word>(clsName));
}


Foam::wordList Foam::IOobjectList::names
(
    const char* clsName,
    const bool syncPar
) const
{
    // No nullptr check - only called with string literals
    return sortedNames(static_cast<word>(clsName), syncPar);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::wordList Foam::IOobjectList::sortedNames() const
{
    return HashPtrTable<IOobject>::sortedToc();
}


Foam::wordList Foam::IOobjectList::sortedNames(const bool syncPar) const
{
    wordList objNames(HashPtrTable<IOobject>::sortedToc());

    checkNameOrder(objNames, syncPar);
    return objNames;
}


Foam::wordList Foam::IOobjectList::sortedNames(const char* clsName) const
{
    // No nullptr check - only called with string literals
    return sortedNames(static_cast<word>(clsName));
}


Foam::wordList Foam::IOobjectList::sortedNames
(
    const char* clsName,
    const bool syncPar
) const
{
    // No nullptr check - only called with string literals
    return sortedNames(static_cast<word>(clsName), syncPar);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::label Foam::IOobjectList::prune_0()
{
    return
        HashPtrTable<IOobject>::filterKeys
        (
            [](const word& k){ return k.ends_with("_0"); },
            true  // prune
        );
}


Foam::wordList Foam::IOobjectList::allNames() const
{
    wordList objNames(HashPtrTable<IOobject>::toc());

    syncNames(objNames);
    return objNames;
}


void Foam::IOobjectList::checkNames(const bool syncPar) const
{
    if (syncPar && UPstream::is_parallel())
    {
        wordList objNames(HashPtrTable<IOobject>::sortedToc());

        checkNameOrder(objNames, syncPar);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const IOobjectList& list)
{
    os << nl << list.size() << nl << token::BEGIN_LIST << nl;

    forAllConstIters(list, iter)
    {
        os  << iter.key() << token::SPACE
            << iter.val()->headerClassName() << nl;
    }

    os << token::END_LIST;
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
