/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2023 OpenCFD Ltd.
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

#include "faBoundaryMeshEntries.H"
#include "processorFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(faBoundaryMeshEntries);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::faBoundaryMeshEntries::faBoundaryMeshEntries(const IOobject& io)
:
    regIOobject
    (
        IOobject(io, IOobjectOption::NO_REGISTER)
    )
{
    // readIOcontents()
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        // Read as entries
        Istream& is = readStream(typeName);

        is >> *this;
        close();
    }
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::PtrList<Foam::entry>
Foam::faBoundaryMeshEntries::readContents(const IOobject& io)
{
    faBoundaryMeshEntries reader(io);

    return PtrList<entry>(std::move(static_cast<PtrList<entry>&>(reader)));
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Extract optional entry from dictionaries and return as a list
template<class T>
static inline List<T> extract
(
    const word& key,
    const UPtrList<entry>& entries,
    const T& initValue
)
{
    List<T> result(entries.size(), initValue);

    forAll(entries, i)
    {
        const dictionary& dict = entries[i].dict();
        dict.readIfPresent(key, result[i]);
    }

    return result;
}

}  // End namespace


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::faBoundaryMeshEntries::removeProcPatches(PtrList<entry>& entries)
{
    // Truncate at the first processor patch entry

    label nNonProcessor = entries.size();

    forAll(entries, patchi)
    {
        const dictionary& dict = entries[patchi].dict();

        const word pType = dict.get<word>("type");
        if (pType == processorFaPatch::typeName)
        {
            nNonProcessor = patchi;
            break;
        }
    }

    entries.resize(nNonProcessor);
}


bool Foam::faBoundaryMeshEntries::writeEntries
(
    Ostream& os,
    const UPtrList<entry>& entries
)
{
    os  << entries.size();

    if (entries.empty())
    {
        // 0-sized : can write with less vertical space
        os  << token::BEGIN_LIST << token::END_LIST;
    }
    else
    {
        os  << nl << token::BEGIN_LIST << incrIndent << nl;

        forAll(entries, patchi)
        {
            const auto& key = entries[patchi].keyword();
            const auto& dict = entries[patchi].dict();

            dict.writeEntry(key, os);
        }
        os  << decrIndent << token::END_LIST;
    }

    os.check(FUNCTION_NAME);
    return os.good();
}


Foam::wordList Foam::faBoundaryMeshEntries::types
(
    const UPtrList<entry>& entries
)
{
    return extract<word>("type", entries, "patch");
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

void Foam::faBoundaryMeshEntries::removeProcPatches()
{
    removeProcPatches(*this);
}


Foam::wordList Foam::faBoundaryMeshEntries::types() const
{
    return extract<word>("type", *this, "patch");
}


void Foam::faBoundaryMeshEntries::writeEntry(Ostream& os) const
{
    writeEntries(os, *this);
}


void Foam::faBoundaryMeshEntries::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    const PtrList<entry>& entries = *this;

    if (!keyword.empty())
    {
        os.write(keyword);
        os << (entries.empty() ? token::SPACE : token::NL);
    }

    writeEntries(os, entries);

    if (!keyword.empty()) os.endEntry();
}


bool Foam::faBoundaryMeshEntries::writeData(Ostream& os) const
{
    return writeEntries(os, *this);
}


bool Foam::faBoundaryMeshEntries::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    streamOpt.compression(IOstreamOption::UNCOMPRESSED);
    return regIOobject::writeObject(streamOpt, writeOnProc);
}


// ************************************************************************* //
