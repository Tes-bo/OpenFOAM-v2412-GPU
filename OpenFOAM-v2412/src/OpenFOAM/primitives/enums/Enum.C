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

#include "Enum.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EnumType>
Foam::Enum<EnumType>::Enum
(
    std::initializer_list<std::pair<EnumType, const char*>> list
)
:
    keys_(list.size()),
    vals_(list.size())
{
    label i = 0;
    for (const auto& pair : list)
    {
        keys_[i] = pair.second;
        vals_[i] = int(pair.first);
        ++i;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EnumType>
void Foam::Enum<EnumType>::push_back
(
    std::initializer_list<std::pair<EnumType, const char*>> list
)
{
    label i = size();

    keys_.resize(i + list.size());
    vals_.resize(i + list.size());

    for (const auto& pair : list)
    {
        keys_[i] = pair.second;
        vals_[i] = int(pair.first);
        ++i;
    }
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::get(const word& enumName) const
{
    const label idx = keys_.find(enumName);

    if (idx < 0)
    {
        FatalErrorInFunction
            << enumName << " is not in enumeration: " << *this << nl
            << exit(FatalError);
    }

    return EnumType(vals_[idx]);
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::lookup
(
    const word& enumName,
    const EnumType deflt
) const
{
    const label idx = keys_.find(enumName);

    if (idx < 0)
    {
        return deflt;
    }

    return EnumType(vals_[idx]);
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::read(Istream& is) const
{
    const word enumName(is);

    const label idx = keys_.find(enumName);

    if (idx < 0)
    {
        FatalIOErrorInFunction(is)
            << enumName << " is not in enumeration: " << *this << nl
            << exit(FatalIOError);
    }

    return EnumType(vals_[idx]);
}


template<class EnumType>
bool Foam::Enum<EnumType>::read
(
    Istream& is,
    EnumType& val,
    const bool mandatory
) const
{
    const word enumName(is);

    const label idx = keys_.find(enumName);

    if (idx >= 0)
    {
        val = EnumType(vals_[idx]);
        return true;
    }

    if (mandatory)
    {
        FatalIOErrorInFunction(is)
            << enumName << " is not in enumeration: " << *this << nl
            << exit(FatalIOError);
    }

    return false;
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::get
(
    const word& key,
    const dictionary& dict
) const
{
    const word enumName(dict.get<word>(key, keyType::LITERAL));

    const label idx = keys_.find(enumName);

    if (idx < 0)
    {
        FatalIOErrorInFunction(dict)
            << "Lookup:" << key << " enumeration " << enumName
            << " is not in enumeration: " << *this << nl
            << exit(FatalIOError);
    }

    return EnumType(vals_[idx]);
}


template<class EnumType>
EnumType Foam::Enum<EnumType>::getOrDefault
(
    const word& key,
    const dictionary& dict,
    const EnumType deflt,
    const bool warnOnly
) const
{
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);

    if (eptr)
    {
        const word enumName(eptr->get<word>());

        const label idx = keys_.find(enumName);

        if (idx >= 0)
        {
            return EnumType(vals_[idx]);
        }

        // Found the dictionary entry, but the name not in enumeration

        if (warnOnly)
        {
            IOWarningInFunction(dict)
                << "Lookup:" << key << " enumeration " << enumName
                << " is not in enumeration: " << *this << nl
                << "using default " << get(deflt)
                << " (value " << int(deflt) << ')' << endl;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Lookup:" << key << " enumeration " << enumName
                << " is not in enumeration: " << *this << nl
                << exit(FatalIOError);
        }
    }

    return deflt;
}


template<class EnumType>
bool Foam::Enum<EnumType>::readEntry
(
    const word& key,
    const dictionary& dict,
    EnumType& val,
    const bool mandatory,
    const bool warnOnly
) const
{
    const entry* eptr = dict.findEntry(key, keyType::LITERAL);

    if (eptr)
    {
        const word enumName(eptr->get<word>());

        const label idx = keys_.find(enumName);

        if (idx >= 0)
        {
            val = EnumType(vals_[idx]);
            return true;
        }

        // Found the dictionary entry, but the name not in enumeration

        if (warnOnly)
        {
            IOWarningInFunction(dict)
                << "Lookup:" << key << " enumeration " << enumName
                << " is not in enumeration: " << *this << nl
                << "leaving value unchanged"
                << " (value " << int(val) << ')' << endl;
        }
        else
        {
            FatalIOErrorInFunction(dict)
                << "Lookup:" << key << " enumeration " << enumName
                << " is not in enumeration: " << *this << nl
                << exit(FatalIOError);
        }
    }
    else if (mandatory)
    {
        FatalIOErrorInFunction(dict)
            << "Lookup:" << key
            << " not found in dictionary " << dict.name() << nl
            << exit(FatalIOError);
    }

    return false;
}


// ************************************************************************* //
