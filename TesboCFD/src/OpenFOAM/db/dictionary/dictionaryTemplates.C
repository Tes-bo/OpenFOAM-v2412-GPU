/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
#include "primitiveEntry.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T>
void Foam::dictionary::reportDefault
(
    const word& keyword,
    const T& deflt,
    const bool added
) const
{
    if (writeOptionalEntries > 1)
    {
        FatalIOError(dictionary::executableName(), *this)
            << "No optional entry: " << keyword
            << " Default: " << deflt << nl
            << exit(FatalIOError);
    }

    OSstream& os = InfoErr.stream(reportingOutput.get());

    // Tag with "-- " prefix to make the message stand out
    os  << "-- Executable: "
        << dictionary::executableName()
        << " Dictionary: ";

    // Double-quote dictionary and entry for more reliably parsing,
    // especially if the keyword contains regular expressions.

    if (this->isNullDict())
    {
        // Output as "", but could have "(null)" etc
        os << token::DQUOTE << token::DQUOTE;
    }
    else
    {
        os.writeQuoted(this->relativeName(), true);
    }

    os  << " Entry: ";
    os.writeQuoted(keyword, true);
    os  << " Default: " << deflt;

    if (added)
    {
        os  << " Added: true";
    }
    os  << nl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Compare>
Foam::wordList Foam::dictionary::sortedToc(const Compare& comp) const
{
    return hashedEntries_.sortedToc(comp);
}


template<class T>
Foam::entry* Foam::dictionary::add(const keyType& k, const T& v, bool overwrite)
{
    return add(new primitiveEntry(k, v), overwrite);
}


template<class T>
Foam::entry* Foam::dictionary::set(const keyType& k, const T& v)
{
    return set(new primitiveEntry(k, v));
}


template<class T>
T Foam::dictionary::get
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    T val;
    readEntry<T>(keyword, val, matchOpt);
    return val;
}


template<class T, class Predicate>
T Foam::dictionary::getCheck
(
    const word& keyword,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    T val;
    readCheck<T, Predicate>(keyword, val, pred, matchOpt);
    return val;
}


template<class T>
T Foam::dictionary::getCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    enum keyType::option matchOpt
) const
{
    T val;
    readCompat<T>(keyword, compat, val, matchOpt);
    return val;
}


template<class T>
T Foam::dictionary::getOrDefault
(
    const word& keyword,
    const T& deflt,
    enum keyType::option matchOpt
) const
{
    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        T val;

        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        reportDefault(keyword, deflt);
    }

    return deflt;
}


template<class T>
T Foam::dictionary::getOrAdd
(
    const word& keyword,
    const T& deflt,
    enum keyType::option matchOpt
)
{
    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        T val;

        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        reportDefault(keyword, deflt, true);  // Added
    }

    add(new primitiveEntry(keyword, deflt));
    return deflt;
}


template<class T, class Predicate>
T Foam::dictionary::getCheckOrDefault
(
    const word& keyword,
    const T& deflt,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    #ifdef FULLDEBUG
    if (!pred(deflt))
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' with invalid default in dictionary "
            << name()
            << exit(FatalIOError);
    }
    #endif

    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        T val;

        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return val;
    }
    else if (writeOptionalEntries)
    {
        reportDefault(keyword, deflt);
    }

    return deflt;
}


template<class T, class Predicate>
T Foam::dictionary::getCheckOrAdd
(
    const word& keyword,
    const T& deflt,
    const Predicate& pred,
    enum keyType::option matchOpt
)
{
    #ifdef FULLDEBUG
    if (!pred(deflt))
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' with invalid default in dictionary "
            << name()
            << exit(FatalIOError);
    }
    #endif

    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        T val;

        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return val;
    }
    else if (writeOptionalEntries)
    {
        reportDefault(keyword, deflt, true);  // Added
    }

    add(new primitiveEntry(keyword, deflt));
    return deflt;
}


template<class T>
bool Foam::dictionary::readEntry
(
    const word& keyword,
    T& val,
    enum keyType::option matchOpt,
    IOobjectOption::readOption readOpt
) const
{
    if (readOpt == IOobjectOption::NO_READ)
    {
        return false;
    }

    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }
    else if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T, class Predicate>
bool Foam::dictionary::readCheck
(
    const word& keyword,
    T& val,
    const Predicate& pred,
    enum keyType::option matchOpt,
    IOobjectOption::readOption readOpt
) const
{
    if (readOpt == IOobjectOption::NO_READ)
    {
        return false;
    }

    const entry* eptr = csearch(keyword, matchOpt).ptr();

    if (eptr)
    {
        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        if (!pred(val))
        {
            raiseBadInput(is, keyword);
        }

        return true;
    }
    else if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T>
bool Foam::dictionary::readCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    T& val,
    enum keyType::option matchOpt,
    IOobjectOption::readOption readOpt
) const
{
    if (readOpt == IOobjectOption::NO_READ)
    {
        return false;
    }

    const entry* eptr = csearchCompat(keyword, compat, matchOpt).ptr();

    if (eptr)
    {
        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return true;
    }
    else if (IOobjectOption::isReadRequired(readOpt))
    {
        FatalIOErrorInFunction(*this)
            << "Entry '" << keyword << "' not found in dictionary "
            << name() << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class T>
bool Foam::dictionary::readIfPresent
(
    const word& keyword,
    T& val,
    enum keyType::option matchOpt
) const
{
    // Reading is optional
    return readEntry<T>
    (
        keyword,
        val,
        matchOpt,
        IOobjectOption::READ_IF_PRESENT
    );
}


template<class T, class Predicate>
bool Foam::dictionary::readCheckIfPresent
(
    const word& keyword,
    T& val,
    const Predicate& pred,
    enum keyType::option matchOpt
) const
{
    // Reading is optional
    return readCheck<T, Predicate>
    (
        keyword,
        val,
        pred,
        matchOpt,
        IOobjectOption::READ_IF_PRESENT
    );
}


template<class T>
T Foam::dictionary::getOrDefaultCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    const T& deflt,
    enum keyType::option matchOpt
) const
{
    const entry* eptr = csearchCompat(keyword, compat, matchOpt).ptr();

    if (eptr)
    {
        T val;

        ITstream& is = eptr->stream();
        is >> val;

        checkITstream(is, keyword);

        return val;
    }
    else if (writeOptionalEntries)
    {
        reportDefault(keyword, deflt);
    }

    return deflt;
}


template<class T>
bool Foam::dictionary::readIfPresentCompat
(
    const word& keyword,
    std::initializer_list<std::pair<const char*,int>> compat,
    T& val,
    enum keyType::option matchOpt
) const
{
    // Reading is optional
    return readCompat<T>
    (
        keyword,
        compat,
        val,
        matchOpt,
        IOobjectOption::READ_IF_PRESENT
    );
}


// ************************************************************************* //
