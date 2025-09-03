/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "argList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::argList::readList(ITstream& is, List<T>& list)
{
    if (is.size() == 1)
    {
        // Single token - treat like List with one entry
        list.resize(1);
        is >> list.front();
    }
    else
    {
        is >> list;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::word& Foam::argList::executable() const noexcept
{
    return executable_;
}


inline const Foam::string& Foam::argList::commandLine() const noexcept
{
    return commandLine_;
}


inline const Foam::fileName& Foam::argList::rootPath() const noexcept
{
    return rootPath_;
}


inline const Foam::fileName& Foam::argList::caseName() const noexcept
{
    return case_;
}


inline const Foam::fileName& Foam::argList::globalCaseName() const noexcept
{
    return globalCase_;
}


inline Foam::fileName Foam::argList::path() const
{
    return rootPath()/caseName();
}


inline Foam::fileName Foam::argList::globalPath() const
{
    return rootPath()/globalCaseName();
}


inline Foam::fileName Foam::argList::relativePath
(
    const fileName& input,
    const bool caseTag
) const
{
    return input.relative(globalPath(), caseTag);
}


inline const Foam::ParRunControl&
Foam::argList::runControl() const noexcept
{
    return runControl_;
}


inline bool Foam::argList::distributed() const noexcept
{
    return runControl_.distributed();
}


inline int Foam::argList::dryRun() const noexcept
{
    return runControl_.dryRun();
}


inline int Foam::argList::dryRun(const int level) noexcept
{
    return runControl_.dryRun(level);
}


inline int Foam::argList::verbose() const noexcept
{
    return runControl_.verbose();
}


inline int Foam::argList::verbose(const int level) noexcept
{
    return runControl_.verbose(level);
}


inline Foam::dlLibraryTable& Foam::argList::libs() const noexcept
{
    return libs_;
}


inline Foam::label Foam::argList::size() const noexcept
{
    return args_.size();
}


inline const Foam::stringList& Foam::argList::args() const noexcept
{
    return args_;
}


inline Foam::stringList& Foam::argList::args() noexcept
{
    return args_;
}


inline const Foam::HashTable<Foam::string>&
Foam::argList::options() const noexcept
{
    return options_;
}


inline Foam::HashTable<Foam::string>&
Foam::argList::options() noexcept
{
    return options_;
}


inline bool Foam::argList::found(const word& optName) const
{
    return options_.contains(optName);
}


inline Foam::ITstream Foam::argList::lookup(const word& optName) const
{
    return ITstream(options_[optName]);
}


// * * * * * * * * * * * * Template Specializations  * * * * * * * * * * * * //

namespace Foam
{
    template<> inline int32_t argList::get<int32_t>(const label index) const
    {
        return Foam::readInt32(args_[index]);
    }

    template<> inline int64_t argList::get<int64_t>(const label index) const
    {
        return Foam::readInt64(args_[index]);
    }

    template<> inline float argList::get<float>(const label index) const
    {
        return Foam::readFloat(args_[index]);
    }

    template<> inline double argList::get<double>(const label index) const
    {
        return Foam::readDouble(args_[index]);
    }


    template<> inline int32_t argList::get<int32_t>(const word& optName) const
    {
        return Foam::readInt32(options_[optName]);
    }

    template<> inline int64_t argList::get<int64_t>(const word& optName) const
    {
        return Foam::readInt64(options_[optName]);
    }

    template<> inline float argList::get<float>(const word& optName) const
    {
        return Foam::readFloat(options_[optName]);
    }

    template<> inline double argList::get<double>(const word& optName) const
    {
        return Foam::readDouble(options_[optName]);
    }


    template<>
    inline string argList::get<Foam::string>(const label index) const
    {
        return args_[index];
    }

    template<>
    inline word argList::get<Foam::word>(const label index) const
    {
        return args_[index];
    }

    template<>
    inline fileName argList::get<Foam::fileName>(const label index) const
    {
        return fileName::validate(args_[index]);
    }


    template<>
    inline string argList::get<Foam::string>(const word& optName) const
    {
        return options_[optName];
    }

    template<>
    inline word argList::get<Foam::word>(const word& optName) const
    {
        return options_[optName];
    }

    template<>
    inline fileName argList::get<Foam::fileName>(const word& optName) const
    {
        return fileName::validate(options_[optName]);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline T Foam::argList::get(const label index) const
{
    ITstream is(args_[index]);

    T val;
    is >> val;

    checkITstream(is, index);

    return val;
}


template<class T>
inline T Foam::argList::get(const word& optName) const
{
    ITstream is(options_[optName]);

    T val;
    is >> val;

    checkITstream(is, optName);

    return val;
}


template<class T>
inline T Foam::argList::getOrDefault
(
    const word& optName,
    const T& deflt
) const
{
    if (options_.contains(optName))
    {
        return get<T>(optName);
    }

    return deflt;
}


template<class T>
inline bool Foam::argList::readIfPresent
(
    const word& optName,
    T& val
) const
{
    if (options_.contains(optName))
    {
        val = get<T>(optName);
        return true;
    }

    return false;
}


template<class T>
inline bool Foam::argList::readIfPresent
(
    const word& optName,
    T& val,
    const T& deflt
) const
{
    if (readIfPresent<T>(optName, val))
    {
        return true;
    }

    val = deflt;
    return false;
}


template<class T>
inline Foam::List<T> Foam::argList::getList(const label index) const
{
    ITstream is(args_[index]);

    List<T> list;
    readList(is, list);

    checkITstream(is, index);

    return list;
}


template<class T>
inline Foam::List<T> Foam::argList::getList
(
    const word& optName,
    bool mandatory
) const
{
    List<T> list;

    if (mandatory || options_.contains(optName))
    {
        ITstream is(options_[optName]);

        readList(is, list);

        checkITstream(is, optName);
    }

    return list;
}


template<class T>
inline bool Foam::argList::readListIfPresent
(
    const word& optName,
    List<T>& list
) const
{
    if (options_.contains(optName))
    {
        ITstream is(options_[optName]);

        readList(is, list);

        checkITstream(is, optName);

        return true;
    }

    return false;
}


template<class T, class Predicate>
inline bool Foam::argList::readCheck
(
    const word& optName,
    T& val,
    const Predicate& pred,
    bool mandatory
) const
{
    if (readIfPresent<T>(optName, val))
    {
        if (!pred(val))
        {
            raiseBadInput(optName);
        }

        return true;
    }
    else if (mandatory)
    {
        FatalError(executable())
            << "Option -" << optName << " not specified" << nl
            << exit(FatalError);
    }

    return false;
}


template<class T, class Predicate>
inline bool Foam::argList::readCheckIfPresent
(
    const word& optName,
    T& val,
    const Predicate& pred
) const
{
    return readCheck<T>(optName, val, pred, false);
}


template<class T, class Predicate>
inline T Foam::argList::getCheck
(
    const word& optName,
    const Predicate& pred
) const
{
    T val;
    readCheck<T>(optName, val, pred, true);
    return val;
}


template<class T, class Predicate>
inline T Foam::argList::getCheckOrDefault
(
    const word& optName,
    const T& deflt,
    const Predicate& pred
) const
{
    // Could predicate check default as well (for FULLDEBUG)

    T val;
    if (readCheck<T>(optName, val, pred, false))
    {
        return val;
    }

    return deflt;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline const Foam::string& Foam::argList::operator[](const label index) const
{
    return args_[index];
}


inline const Foam::string& Foam::argList::operator[](const word& optName) const
{
    return options_[optName];
}


// ************************************************************************* //
