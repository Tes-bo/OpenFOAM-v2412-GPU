/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
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

#include "dictionary.H"
#include "error.H"
#include "JobInfo.H"
#include "primitiveEntry.H"
#include "dictionaryEntry.H"
#include "regExp.H"
#include "OSHA1stream.H"
#include "OSstream.H"
#include "argList.H"
#include "registerSwitch.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(dictionary, 0);
}

Foam::refPtr<Foam::OSstream> Foam::dictionary::reportingOutput(nullptr);

const Foam::dictionary Foam::dictionary::null;

int Foam::dictionary::writeOptionalEntries
(
    Foam::debug::infoSwitch("writeOptionalEntries", 0)
);


registerInfoSwitch
(
    "writeOptionalEntries",
    int,
    Foam::dictionary::writeOptionalEntries
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::word Foam::dictionary::executableName()
{
    return argList::envExecutable();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionary::dictionary()
:
    name_(),
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary(const fileName& name)
:
    name_(name),
    parent_(dictionary::null)
{}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    const dictionary& dict
)
:
    parent_type(dict, *this),
    name_(dict.name()),
    parent_(parentDict)
{
    for (entry& e : *this)
    {
        hashedEntries_.insert(e.keyword(), &e);

        if (e.keyword().isPattern())
        {
            patterns_.push_front(&e);
            regexps_.push_front(autoPtr<regExp>::New(e.keyword()));
        }
    }
}


Foam::dictionary::dictionary
(
    const dictionary& dict
)
:
    parent_type(dict, *this),
    name_(dict.name()),
    parent_(dictionary::null)
{
    for (entry& e : *this)
    {
        hashedEntries_.insert(e.keyword(), &e);

        if (e.keyword().isPattern())
        {
            patterns_.push_front(&e);
            regexps_.push_front(autoPtr<regExp>::New(e.keyword()));
        }
    }
}


Foam::dictionary::dictionary(const dictionary* dict)
:
    name_(),
    parent_(dictionary::null)
{
    if (dict)
    {
        operator=(*dict);
    }
}


Foam::dictionary::dictionary
(
    const dictionary& parentDict,
    dictionary&& dict
)
:
    name_(),
    parent_(parentDict)
{
    transfer(dict);
    name() = fileName::concat(parentDict.name(), name(), '/');
}


Foam::dictionary::dictionary
(
    dictionary&& dict
)
:
    name_(),
    parent_(dictionary::null)
{
    transfer(dict);
}


Foam::autoPtr<Foam::dictionary> Foam::dictionary::clone() const
{
    return autoPtr<dictionary>::New(*this);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dictionary::~dictionary()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::fileName Foam::dictionary::relativeName(const bool caseTag) const
{
    return argList::envRelativePath(name(), caseTag);
}


const Foam::dictionary& Foam::dictionary::topDict() const
{
    const dictionary& p = parent();

    if (&p != this && !p.name().empty())
    {
        return p.topDict();
    }

    return *this;
}


Foam::label Foam::dictionary::startLineNumber() const
{
    return
    (
        IDLList<entry>::empty()
      ? -1
      : IDLList<entry>::front()->startLineNumber()
    );
}


Foam::label Foam::dictionary::endLineNumber() const
{
    return
    (
        IDLList<entry>::empty()
      ? -1
      : IDLList<entry>::back()->endLineNumber()
    );
}


Foam::SHA1Digest Foam::dictionary::digest() const
{
    OSHA1stream os;

    // Process entries
    for (const entry& e : *this)
    {
        os << e;
    }

    return os.digest();
}


Foam::tokenList Foam::dictionary::tokens() const
{
    // Serialize dictionary entries into a string
    OCharStream os;

    // Process entries
    for (const entry& e : *this)
    {
        os << e;
    }

    // String re-parsed as a list of tokens
    return ITstream::parse(os.view());
}


void Foam::dictionary::checkITstream
(
    const ITstream& is,
    const word& keyword
) const
{
    const label remaining = (is.size() ? is.nRemainingTokens() : -100);

    if (!remaining)
    {
        return;
    }

    // Similar to SafeFatalIOError
    if (JobInfo::constructed)
    {
        OSstream& err =
            FatalIOError
            (
                "",                 // functionName
                "",                 // sourceFileName
                0,                  // sourceFileLineNumber
                relativeName(),     // ioFileName == dictionary name
                is.lineNumber()     // ioStartLineNumber
            );

        if (remaining > 0)
        {
            err
                << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl
                << "    ";
            is.writeList(err, 0);  // <- flatOutput
        }
        else
        {
            err << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl;
        }

        err << exit(FatalIOError);
    }
    else
    {
        // Not yet constructed

        std::cerr
            << nl
            << "--> FOAM FATAL IO ERROR:" << nl;

        if (remaining > 0)
        {
            std::cerr
                << "Entry '" << keyword << "' has "
                << remaining << " excess tokens in stream" << nl << nl;
        }
        else
        {
            std::cerr
                << "Entry '" << keyword
                << "' had no tokens in stream" << nl << nl;
        }

        std::cerr
            // ioFileName == dictionary name
            << "file: " << relativeName()
            << " at line " << is.lineNumber() << '.' << nl
            << std::endl;

        std::exit(1);
    }
}


void Foam::dictionary::raiseBadInput
(
    const ITstream& is,
    const word& keyword
) const
{
    // Can use FatalIOError instead of SafeFatalIOError
    // since predicate checks are not used at the earliest stages
    FatalIOError
    (
        "",                 // functionName
        "",                 // sourceFileName
        0,                  // sourceFileLineNumber
        relativeName(),     // ioFileName == dictionary name
        is.lineNumber(),    // ioStartLineNumber
        -1                  // ioEndLineNumber
    )
        << "Entry '" << keyword << "' with invalid input" << nl
        << exit(FatalIOError);
}


const Foam::entry& Foam::dictionary::lookupEntry
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    const entry* eptr = findEntry(keyword, matchOpt);

    if (!eptr)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    return *eptr;
}


Foam::ITstream& Foam::dictionary::lookup
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return lookupEntry(keyword, matchOpt).stream();
}


bool Foam::dictionary::substituteKeyword
(
    const word& keyword,
    bool mergeEntry
)
{
    if (keyword.size() < 2)
    {
        return false;
    }

    // Drop leading '$' to get the var-name, already validated as word.
    const word varName(keyword.substr(1), false);

    // Lookup the variable name in the given dictionary
    const const_searcher finder(csearch(varName, keyType::REGEX_RECURSIVE));

    // If defined insert its entries into this dictionary
    if (finder.good())
    {
        for (const entry& e : finder.dict())
        {
            add(e, mergeEntry);
        }

        return true;
    }

    return false;
}


bool Foam::dictionary::substituteScopedKeyword
(
    const word& keyword,
    bool mergeEntry
)
{
    if (keyword.size() < 2)
    {
        return false;
    }

    // Drop leading '$' to get the var-name, already validated as word.
    const word varName(keyword.substr(1), false);

    // Lookup the variable name in the given dictionary
    const auto finder(csearchScoped(varName, keyType::REGEX_RECURSIVE));

    // If defined insert its entries into this dictionary
    if (finder.good())
    {
        for (const entry& e : finder.dict())
        {
            add(e, mergeEntry);
        }

        return true;
    }

    return false;
}


const Foam::dictionary& Foam::dictionary::subDict
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    if (!finder.good())
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    return finder.dict();
}


Foam::dictionary& Foam::dictionary::subDict
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    searcher finder(search(keyword, matchOpt));

    if (!finder.good())
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    return finder.dict();
}


Foam::dictionary& Foam::dictionary::subDictOrAdd
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    searcher finder(search(keyword, matchOpt));

    dictionary* dictPtr = finder.dictPtr();

    if (dictPtr)
    {
        // Found and a sub-dictionary
        return *dictPtr;
    }

    if (finder.good())
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword
            << "' is not a sub-dictionary in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    dictPtr = this->set(keyword, dictionary())->dictPtr();

    if (!dictPtr)
    {
        FatalIOErrorInFunction(*this)
            << "Failed to insert sub-dictionary '" << keyword
            << "' in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    return *dictPtr;
}


Foam::dictionary Foam::dictionary::subOrEmptyDict
(
    const word& keyword,
    enum keyType::option matchOpt,
    const bool mandatory
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    const dictionary* dictPtr = finder.dictPtr();

    if (dictPtr)
    {
        // Found and a sub-dictionary
        return *dictPtr;
    }

    if (mandatory)
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword
            << "' is not a sub-dictionary in dictionary "
            << relativeName() << nl
            << exit(FatalIOError);
    }

    if (finder.good())
    {
        IOWarningInFunction(*this)
            << "Entry '" << keyword
            << "' found but not a sub-dictionary in dictionary "
            << relativeName() << endl;
    }

    // The move constructor properly qualifies the dictionary name
    return dictionary(*this, dictionary(fileName(keyword)));
}


const Foam::dictionary& Foam::dictionary::optionalSubDict
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    const const_searcher finder(csearch(keyword, matchOpt));

    const dictionary* dictPtr = finder.dictPtr();

    if (dictPtr)
    {
        // Found and a sub-dictionary
        return *dictPtr;
    }

    if (finder.good())
    {
        IOWarningInFunction(*this)
            << "Entry '" << keyword
            << "' found but not a sub-dictionary in dictionary "
            << relativeName() << endl;
    }

    return *this;
}


Foam::wordList Foam::dictionary::toc() const
{
    wordList list(size());

    label n = 0;
    for (const entry& e : *this)
    {
        list[n++] = e.keyword();
    }

    return list;
}


Foam::wordList Foam::dictionary::sortedToc() const
{
    return hashedEntries_.sortedToc();
}


Foam::List<Foam::keyType> Foam::dictionary::keys(bool patterns) const
{
    List<keyType> list(size());

    label n = 0;
    for (const entry& e : *this)
    {
        if (e.keyword().isPattern() ? patterns : !patterns)
        {
            list[n++] = e.keyword();
        }
    }
    list.resize(n);

    return list;
}


Foam::entry* Foam::dictionary::add(entry* entryPtr, bool mergeEntry)
{
    if (!entryPtr)
    {
        return nullptr;
    }

    auto iter = hashedEntries_.find(entryPtr->keyword());

    if (mergeEntry && iter.good())
    {
        // Merge dictionary with dictionary
        if (iter()->isDict() && entryPtr->isDict())
        {
            iter()->dict().merge(entryPtr->dict());

            delete entryPtr;
            return iter();   // pointer to existing dictionary
        }


        // Replace existing dictionary with entry or vice versa
        parent_type::replace(iter(), entryPtr);
        delete iter();
        hashedEntries_.erase(iter);

        if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
        {
            entryPtr->name() =
                fileName::concat(name(), entryPtr->keyword(), '/');

            if (entryPtr->keyword().isPattern())
            {
                patterns_.push_front(entryPtr);
                regexps_.push_front(autoPtr<regExp>::New(entryPtr->keyword()));
            }

            return entryPtr;  // now an entry in the dictionary
        }


        IOWarningInFunction(*this)
            << "Problem replacing entry "<< entryPtr->keyword()
            << " in dictionary " << relativeName() << endl;

        parent_type::remove(entryPtr);

        delete entryPtr;
        return nullptr;
    }


    if (hashedEntries_.insert(entryPtr->keyword(), entryPtr))
    {
        entryPtr->name() =
            fileName::concat(name(), entryPtr->keyword(), '/');

        parent_type::push_back(entryPtr);

        if (entryPtr->keyword().isPattern())
        {
            patterns_.push_front(entryPtr);
            regexps_.push_front(autoPtr<regExp>::New(entryPtr->keyword()));
        }

        return entryPtr;  // now an entry in the dictionary
    }


    IOWarningInFunction(*this)
        << "Attempt to add entry " << entryPtr->keyword()
        << " which already exists in dictionary "
        << relativeName() << endl;

    delete entryPtr;
    return nullptr;
}


Foam::entry* Foam::dictionary::add(const entry& e, bool mergeEntry)
{
    return add(e.clone(*this).ptr(), mergeEntry);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const word& v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const Foam::string& v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const label v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const scalar v,
    bool overwrite
)
{
    return add(new primitiveEntry(k, token(v)), overwrite);
}


Foam::entry* Foam::dictionary::add
(
    const keyType& k,
    const dictionary& v,
    bool mergeEntry
)
{
    return add(new dictionaryEntry(k, *this, v), mergeEntry);
}


Foam::entry* Foam::dictionary::set(entry* entryPtr)
{
    if (!entryPtr)
    {
        return nullptr;
    }

    // Find non-recursive with patterns
    searcher finder(search(entryPtr->keyword(), keyType::REGEX));

    dictionary* dictPtr = finder.dictPtr();

    // Clear dictionary so merge acts like overwrite
    if (dictPtr)
    {
        dictPtr->clear();
    }

    return add(entryPtr, true);
}


Foam::entry* Foam::dictionary::set(const entry& e)
{
    return set(e.clone(*this).ptr());
}


Foam::entry* Foam::dictionary::set(const keyType& k, const dictionary& v)
{
    return set(new dictionaryEntry(k, *this, v));
}


bool Foam::dictionary::merge(const dictionary& dict)
{
    if (this == &dict)
    {
        FatalIOErrorInFunction(*this)
            << "Attempted merge to self, for dictionary "
            << relativeName() << nl
            << abort(FatalIOError);
    }

    bool changed = false;

    for (const entry& e : dict)
    {
        auto fnd = hashedEntries_.find(e.keyword());

        if (fnd.good())
        {
            // Recursively merge sub-dictionaries
            // TODO: merge without copying
            if (fnd()->isDict() && e.isDict())
            {
                if (fnd()->dict().merge(e.dict()))
                {
                    changed = true;
                }
            }
            else
            {
                add(e.clone(*this).ptr(), true);
                changed = true;
            }
        }
        else
        {
            // Not found - just add
            add(e.clone(*this).ptr());
            changed = true;
        }
    }

    return changed;
}


void Foam::dictionary::clear()
{
    parent_type::clear();
    hashedEntries_.clear();
    patterns_.clear();
    regexps_.clear();
}


void Foam::dictionary::transfer(dictionary& dict)
{
    // Changing parents probably doesn't make much sense,
    // but what about the names?
    name() = dict.name();

    parent_type::transfer(dict);
    hashedEntries_.transfer(dict.hashedEntries_);
    patterns_.transfer(dict.patterns_);
    regexps_.transfer(dict.regexps_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::dictionary::operator=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        return;  // Self-assignment is a no-op
    }

    name() = rhs.name();
    clear();

    // Create clones of the entries in the given dictionary
    // resetting the parentDict to this dictionary

    for (const entry& e : rhs)
    {
        add(e.clone(*this).ptr());
    }
}


void Foam::dictionary::operator+=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "Attempted addition to self, for dictionary "
            << relativeName() << nl
            << abort(FatalIOError);
    }

    for (const entry& e : rhs)
    {
        add(e.clone(*this).ptr());
    }
}


void Foam::dictionary::operator|=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "Attempted |= merging to self, for dictionary "
            << relativeName() << nl
            << abort(FatalIOError);
    }

    for (const entry& e : rhs)
    {
        if (!found(e.keyword()))
        {
            add(e.clone(*this).ptr());
        }
    }
}


void Foam::dictionary::operator<<=(const dictionary& rhs)
{
    if (this == &rhs)
    {
        FatalIOErrorInFunction(*this)
            << "Attempted addition to self, for dictionary "
            << relativeName() << nl
            << abort(FatalIOError);
    }

    for (const entry& e : rhs)
    {
        set(e.clone(*this).ptr());
    }
}


/* * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * */

Foam::dictionary Foam::operator+
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary result(dict1);
    result += dict2;
    return result;
}


Foam::dictionary Foam::operator|
(
    const dictionary& dict1,
    const dictionary& dict2
)
{
    dictionary result(dict1);
    result |= dict2;
    return result;
}


// ************************************************************************* //
