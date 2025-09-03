/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "colourTable.H"
#include "etcFiles.H"
#include "IFstream.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::HashPtrTable<Foam::colourTable> Foam::colourTable::tables_;



const Foam::Enum
<
    Foam::colourTable::predefinedType
>
Foam::colourTable::predefinedNames
({
    { predefinedType::COOL_WARM, "coolToWarm" },
    { predefinedType::COLD_HOT, "coldAndHot" },
    { predefinedType::FIRE, "fire" },
    { predefinedType::RAINBOW, "rainbow" },
    { predefinedType::GREYSCALE, "greyscale" },
    { predefinedType::XRAY, "xray" },
});


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::colourTable::constructTables()
{
    if (tables_.size())
    {
        FatalErrorInFunction
            << "attempt to re-construct colourTables when they already exist"
            << exit(FatalError);
    }

    IFstream is(findEtcFile("colourTables", true));  // Mandatory file

    HashPtrTable<colourTable> newEntries(is);
    tables_.swap(newEntries);

    Info<< "loaded " << tables_.sortedToc()
        << " from etc/colourTable" << endl;

    Info<< "== " << tables_ << nl;
}


const Foam::HashPtrTable<Foam::colourTable>& Foam::colourTable::tables()
{
    if (tables_.empty())
    {
        constructTables();
    }

    return tables_;
}


const Foam::colourTable* Foam::colourTable::ptr(const word& tableName)
{
    if (tables_.empty())
    {
        constructTables();
    }

    const auto iter = tables_.cfind(tableName);

    if (iter.good())
    {
        const colourTable* p = iter.val();
        return p;
    }

    return nullptr;
}


const Foam::colourTable* Foam::colourTable::ptr(const predefinedType tbl)
{
    return ptr(predefinedNames[tbl]);
}


const Foam::colourTable& Foam::colourTable::ref(const word& tableName)
{
    const colourTable* p = ptr(tableName);

    if (!p)
    {
        FatalErrorInFunction
            << "No such colourTable: " << tableName
            << exit(FatalError);
    }

    return *p;
}


const Foam::colourTable& Foam::colourTable::ref(const predefinedType tbl)
{
    return ref(predefinedNames[tbl]);
}


// ************************************************************************* //
