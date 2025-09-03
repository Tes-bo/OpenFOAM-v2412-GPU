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

#include "IOobjectList.H"
#include "vtkDataArraySelection.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class StringType>
Foam::label Foam::foamPvCore::addToArray
(
    vtkDataArraySelection *select,
    const std::string& prefix,
    const UList<StringType>& names
)
{
    if (prefix.empty())
    {
        for (const auto& name : names)
        {
            select->AddArray(name.c_str());
        }
    }
    else
    {
        for (const auto& name : names)
        {
            select->AddArray((prefix + name).c_str());
        }
    }

    return names.size();
}


template<class StringType>
Foam::label Foam::foamPvCore::addToArray
(
    vtkDataArraySelection *select,
    const UList<StringType>& names,
    const std::string& suffix
)
{
    if (suffix.empty())
    {
        for (const auto& name : names)
        {
            select->AddArray(name.c_str());
        }
    }
    else
    {
        for (const auto& name : names)
        {
            select->AddArray((name + suffix).c_str());
        }
    }

    return names.size();
}


template<class Type>
Foam::label Foam::foamPvCore::addToSelection
(
    vtkDataArraySelection *select,
    const std::string& prefix,
    const IOobjectList& objects
)
{
    return addToArray(select, prefix, objects.sortedNames(Type::typeName));
}


template<class Type>
Foam::label Foam::foamPvCore::addToSelection
(
    vtkDataArraySelection *select,
    const IOobjectList& objects,
    const std::string& suffix
)
{
    return addToArray(select, objects.sortedNames(Type::typeName), suffix);
}


template<class Type>
Foam::label Foam::foamPvCore::addToSelection
(
    vtkDataArraySelection *select,
    const std::string& prefix,
    const HashTable<wordHashSet>& objects
)
{
    auto iter = objects.cfind(Type::typeName);

    if (iter.good())
    {
        return addToArray(select, prefix, iter.val().sortedToc());
    }

    return 0;
}


template<class Type>
Foam::label Foam::foamPvCore::addToSelection
(
    vtkDataArraySelection *select,
    const HashTable<wordHashSet>& objects,
    const std::string& suffix
)
{
    auto iter = objects.cfind(Type::typeName);

    if (iter.good())
    {
        return addToArray(select, iter.val().sortedToc(), suffix);
    }

    return 0;
}


template<class AnyValue, class AnyHasher>
void Foam::foamPvCore::setSelectedArrayEntries
(
    vtkDataArraySelection* select,
    const HashTable<AnyValue, string, AnyHasher>& enabled
)
{
    const int n = select->GetNumberOfArrays();

    // Disable everything not explicitly enabled
    select->DisableAllArrays();

    // Loop through entries, enabling as required
    for (int i=0; i < n; ++i)
    {
        const char* arrayName = select->GetArrayName(i);
        if (enabled.contains(arrayName))
        {
            select->EnableArray(arrayName);
        }
    }
}


// ************************************************************************* //
