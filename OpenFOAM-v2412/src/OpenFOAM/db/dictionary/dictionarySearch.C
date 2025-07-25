/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "dictionary.H"
#include "dictionaryEntry.H"
#include "stringOps.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{
    // Walk lists of patterns and regexps for an exact match
    // or a regular expression match
    template<class WcIterator, class ReIterator>
    bool findInPatterns
    (
        const bool literal,
        const Foam::word& keyword,
        WcIterator& wcIter,
        ReIterator& reIter
    )
    {
        while (wcIter.good())
        {
            if
            (
                literal
              ? wcIter()->keyword() == keyword
              : reIter()->match(keyword)
            )
            {
                return true;
            }

            ++reIter;
            ++wcIter;
        }

        return false;
    }

} // End anonymous namespace


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearchDotScoped
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    auto scopePos = keyword.find('.');

    if (scopePos == std::string::npos)
    {
        // Normal, non-scoped search
        return csearch(keyword, matchOpt);
    }

    // It is '.' scoped - force non-recusive searching
    matchOpt = keyType::option(matchOpt & ~(keyType::RECURSIVE));

    if (scopePos == 0)
    {
        // Starting with a '.' -> go up for every further '.' found
        ++scopePos;

        const dictionary* dictPtr = this;
        for
        (
            string::const_iterator it = keyword.begin()+1;
            it != keyword.end() && *it == '.';
            ++scopePos, ++it
        )
        {
            // Go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction(*this)
                    << "No parent of current dictionary when searching for "
                    << keyword.substr(1)
                    << exit(FatalIOError);

                return nullptr;
            }
        }

        return dictPtr->csearchDotScoped
        (
            keyword.substr(scopePos),
            matchOpt
        );
    }

    // The first word
    const_searcher finder = csearchDotScoped
    (
        keyword.substr(0, scopePos),
        matchOpt
    );

    // Fall back to finding key with '.' so e.g. if keyword is
    // a.b.c.d it would try
    // a.b, a.b.c, a.b.c.d

    if (!finder.good())
    {
        while (!finder.isDict())
        {
            scopePos = keyword.find('.', scopePos+1);

            // Local entry:
            finder = csearch(keyword.substr(0, scopePos), matchOpt);

            if (scopePos == std::string::npos)
            {
                // Parsed the whole word. Return entry or null.
                return finder;
            }
        }
    }

    if (finder.isDict())
    {
        return finder.dict().csearchDotScoped
        (
            keyword.substr(scopePos),
            matchOpt
        );
    }

    return finder;
}


Foam::dictionary::const_searcher Foam::dictionary::csearchSlashScoped
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{

    // With '/' scoping - recursive is never allowed
    matchOpt = keyType::option(matchOpt & ~(keyType::RECURSIVE));

    const dictionary* dictPtr = this;

    const auto slash = keyword.find('/');

    if (slash == std::string::npos)
    {
        // No slashes:
        // Can use normal (non-scoped) search at the current dictionary level
        return csearch(keyword, matchOpt);
    }
    else if (slash == 0)
    {
        // isAbsolute:
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }
    }

    // Split on '/'
    auto cmpts = stringOps::split<std::string>(keyword, '/');
    auto remaining = cmpts.size();

    for (const auto& cmpt : cmpts)
    {
        --remaining; // Decrement now so we can check (remaining == 0)

        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction(*dictPtr)
                    << "No parent of current dictionary when searching for "
                    << keyword << " at " << cmpt
                    << exit(FatalIOError);
                break;
            }
        }
        else
        {
            // Find entry
            const word key = word::validate(cmpt);

            auto finder = dictPtr->csearch(key, matchOpt);

            if (finder.good())
            {
                if (remaining)
                {
                    // Intermediate must be a dictionary
                    if (finder.isDict())
                    {
                        dictPtr = finder.dictPtr();
                    }
                    else
                    {
                        return const_searcher(dictPtr);
                    }
                }
                else
                {
                    // Last entry - done
                    return finder;
                }
            }
            else
            {
                break;
            }
        }
    }

    // Failed at this dictionary level
    return const_searcher(dictPtr);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::dictionary::const_searcher Foam::dictionary::csearch
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    const_searcher finder(this);

    auto iter = hashedEntries_.cfind(keyword);

    if (iter.good())
    {
        finder.set(iter.val());
        return finder;
    }

    if ((matchOpt & keyType::REGEX) && patterns_.size())
    {
        auto wcLink = patterns_.cbegin();
        auto reLink = regexps_.cbegin();

        // Find in patterns : non-literal matching
        if (findInPatterns(false, keyword, wcLink, reLink))
        {
            finder.set(*wcLink);
            return finder;
        }
    }

    if ((matchOpt & keyType::RECURSIVE) && &parent_ != &dictionary::null)
    {
        return parent_.csearch(keyword, matchOpt);
    }

    return finder;
}


Foam::dictionary::const_searcher Foam::dictionary::search
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearch(keyword, matchOpt);
}


Foam::dictionary::searcher Foam::dictionary::search
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    const_searcher finder = csearch(keyword, matchOpt);

    return static_cast<const searcher&>(finder);
}


Foam::dictionary::const_searcher Foam::dictionary::csearchScoped
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    if (keyword.contains('/'))
    {
        return csearchSlashScoped(keyword, matchOpt);
    }

    if (keyword[0] == ':' || keyword[0] == '^')
    {
        // It is ':' scoped - force non-recursive searching
        matchOpt = keyType::option(matchOpt & ~(keyType::RECURSIVE));

        // Ascend to top-level
        const dictionary* dictPtr = this;
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }

        return dictPtr->csearchDotScoped(keyword.substr(1), matchOpt);
    }

    return csearchDotScoped(keyword, matchOpt);
}


Foam::dictionary::const_searcher Foam::dictionary::searchScoped
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearchScoped(keyword, matchOpt);
}


Foam::dictionary::searcher Foam::dictionary::searchScoped
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    const_searcher finder = csearchScoped(keyword, matchOpt);

    return static_cast<const searcher&>(finder);
}


const Foam::dictionary* Foam::dictionary::cfindScopedDict
(
    const fileName& dictPath
) const
{
    // Or warning
    if (dictPath.empty())
    {
        return nullptr;
    }

    const dictionary* dictPtr = this;
    if (dictPath[0] == '/')
    {
        // isAbsolute:
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = &dictPtr->parent_;
        }
    }

    fileName path(dictPath); // Work on copy
    path.clean();  // Remove unneeded ".."
    auto dictCmpts = stringOps::split(path, '/');  // Split on '/'

    for (const auto& cmpt : dictCmpts)
    {
        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = &dictPtr->parent_;
            }
            else
            {
                FatalIOErrorInFunction(*dictPtr)
                    << "No parent for dictionary while searching "
                    << path
                    << exit(FatalIOError);

                return nullptr;
            }
        }
        else
        {
            // Non-recursive, no patternMatch:
            // do direct lookup, without csearch(cmpt, keyType::LITERAL)

            const word cmptName(cmpt.str(), false);

            auto iter = dictPtr->hashedEntries_.cfind(cmptName);

            if (iter.good())
            {
                const entry *eptr = iter.val();

                if (eptr->isDict())
                {
                    dictPtr = eptr->dictPtr();
                }
                else
                {
                    FatalIOErrorInFunction(*dictPtr)
                        << "Found entry '" << cmptName
                        << "' but not a dictionary, while searching scoped"
                        << nl
                        << "    " << path
                        << exit(FatalIOError);

                    return nullptr;
                }
            }
            else
            {
                return nullptr;
            }
        }
    }

    return dictPtr;
}


const Foam::dictionary* Foam::dictionary::findScopedDict
(
    const fileName& dictPath
) const
{
    return cfindScopedDict(dictPath);
}


Foam::dictionary* Foam::dictionary::findScopedDict
(
    const fileName& dictPath
)
{
    const dictionary* ptr = cfindScopedDict(dictPath);
    return const_cast<dictionary*>(ptr);
}


Foam::dictionary* Foam::dictionary::makeScopedDict(const fileName& dictPath)
{
    // Or warning
    if (dictPath.empty())
    {
        return nullptr;
    }

    dictionary* dictPtr = this;
    if (dictPath[0] == '/')
    {
        // isAbsolute:
        // Ascend to top-level
        while (&dictPtr->parent_ != &dictionary::null)
        {
            dictPtr = const_cast<dictionary*>(&dictPtr->parent_);
        }
    }

    std::string path(dictPath); // Work on a copy
    fileName::clean(path);  // Remove unneeded ".."
    auto dictCmpts = stringOps::split(path, '/'); // Split on '/'

    for (const auto& cmpt : dictCmpts)
    {
        if (cmpt == ".")
        {
            // "." - ignore
        }
        else if (cmpt == "..")
        {
            // ".." - go to parent
            if (&dictPtr->parent_ != &dictionary::null)
            {
                dictPtr = const_cast<dictionary*>(&dictPtr->parent_);
            }
            else
            {
                FatalIOErrorInFunction(*dictPtr)
                    << "No parent for dictionary while searching "
                    << path
                    << exit(FatalIOError);

                return nullptr;
            }
        }
        else
        {
            // Non-recursive, no patternMatch:
            // do direct lookup, without csearch(cmptName, keyType::LITERAL)

            const word cmptName(cmpt.str(), false);

            auto iter = dictPtr->hashedEntries_.find(cmptName);

            if (iter.good())
            {
                entry *eptr = iter.val();

                if (eptr->isDict())
                {
                    dictPtr = eptr->dictPtr();
                }
                else
                {
                    FatalIOErrorInFunction(*dictPtr)
                        << "Cannot create sub-dictionary entry '" << cmptName
                        << "' - a non-dictionary entry is in the way"
                        << nl << "Encountered in scope" << nl
                        << "    " << path
                        << exit(FatalIOError);

                    return nullptr;
                }
            }
            else
            {
                dictionaryEntry *eptr =
                    new dictionaryEntry(cmptName, *dictPtr, dictionary());

                // Add *without* merging, since we just checked that the entry
                // doesn't exist and to ensure that the pointer remains valid.

                if (dictPtr->add(eptr, false))  // NO merge
                {
                    dictPtr = eptr;
                }
                else
                {
                    // Note: a failed add() deletes the eptr passed
                    return nullptr;
                }
            }
        }
    }

    return dictPtr;
}


bool Foam::dictionary::remove(const word& keyword)
{
    auto iter = hashedEntries_.find(keyword);

    if (iter.good())
    {
        // Delete from patterns
        auto wcLink = patterns_.begin();
        auto reLink = regexps_.begin();

        // Find in patterns : literal matching
        if (findInPatterns(true, keyword, wcLink, reLink))
        {
            patterns_.remove(wcLink);
            regexps_.remove(reLink);
        }

        parent_type::remove(iter());
        delete iter();
        hashedEntries_.erase(iter);

        return true;
    }

    return false;
}


bool Foam::dictionary::changeKeyword
(
    const keyType& oldKeyword,
    const keyType& newKeyword,
    bool overwrite
)
{
    // No change
    if (oldKeyword == newKeyword)
    {
        return false;
    }

    // Check that oldKeyword exists and can be changed
    auto iter = hashedEntries_.find(oldKeyword);

    if (!iter.good())
    {
        return false;
    }

    if (iter()->keyword().isPattern())
    {
        FatalIOErrorInFunction(*this)
            << "Old keyword " << oldKeyword << " is a pattern." << nl
            << "Pattern replacement is not supported." << nl
            << exit(FatalIOError);
    }


    auto iter2 = hashedEntries_.find(newKeyword);

    // newKeyword already exists
    if (iter2.good())
    {
        if (overwrite)
        {
            if (iter2()->keyword().isPattern())
            {
                // Delete from patterns
                auto wcLink = patterns_.begin();
                auto reLink = regexps_.begin();

                // Find in patterns : literal matching
                if (findInPatterns(true, iter2()->keyword(), wcLink, reLink))
                {
                    patterns_.remove(wcLink);
                    regexps_.remove(reLink);
                }
            }

            parent_type::replace(iter2(), iter());
            delete iter2();
            hashedEntries_.erase(iter2);
        }
        else
        {
            IOWarningInFunction(*this)
                << "Cannot rename keyword " << oldKeyword
                << " to existing keyword " << newKeyword
                << " in dictionary " << name() << endl;
            return false;
        }
    }

    // Change name and HashTable, but leave DL-List untouched
    iter()->keyword() = newKeyword;
    iter()->name() = fileName::concat(name(), newKeyword, '/');
    hashedEntries_.erase(oldKeyword);
    hashedEntries_.insert(newKeyword, iter());

    if (newKeyword.isPattern())
    {
        patterns_.push_front(iter());
        regexps_.push_front(autoPtr<regExp>::New(newKeyword));
    }

    return true;
}


// ************************************************************************* //
