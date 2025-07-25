/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::dictionary::reportOptional() noexcept
{
    return writeOptionalEntries;
}


inline int Foam::dictionary::reportOptional(const int level) noexcept
{
    int old(writeOptionalEntries);
    writeOptionalEntries = level;
    return old;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fileName& Foam::dictionary::name() const noexcept
{
    return name_;
}


inline Foam::fileName& Foam::dictionary::name() noexcept
{
    return name_;
}


inline Foam::word Foam::dictionary::dictName() const
{
    // With other (non-slash) separator. Eg, with '.'
    // word scopedName(name_.name());
    //
    // const auto i = scopedName.rfind('.');
    // if (i == std::string::npos)
    // {
    //     return scopedName;
    // }
    //
    // return scopedName.substr(i+1);

    // With '/' separator, this is just fileName::name()
    return name_.name();
}


inline bool Foam::dictionary::isNullDict() const noexcept
{
    return (this == &dictionary::null);
}


inline const Foam::dictionary& Foam::dictionary::parent() const noexcept
{
    return parent_;
}


inline const Foam::entry* Foam::dictionary::findEntry
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearch(keyword, matchOpt).ptr();
}


inline Foam::entry* Foam::dictionary::findEntry
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    return const_cast<entry*>(csearch(keyword, matchOpt).ptr());
}


inline bool Foam::dictionary::found
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return static_cast<bool>(findEntry(keyword, matchOpt));
}


inline const Foam::entry* Foam::dictionary::findScoped
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearchScoped(keyword, matchOpt).ptr();
}


inline const Foam::dictionary* Foam::dictionary::findDict
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearch(keyword, matchOpt).dictPtr();
}


inline Foam::dictionary* Foam::dictionary::findDict
(
    const word& keyword,
    enum keyType::option matchOpt
)
{
    return const_cast<dictionary*>(csearch(keyword, matchOpt).dictPtr());
}


inline Foam::ITstream* Foam::dictionary::findStream
(
    const word& keyword,
    enum keyType::option matchOpt
) const
{
    return csearch(keyword, matchOpt).streamPtr();
}


// ************************************************************************* //
