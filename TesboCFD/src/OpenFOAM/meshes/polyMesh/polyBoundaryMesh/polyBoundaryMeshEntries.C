/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "polyBoundaryMeshEntries.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(polyBoundaryMeshEntries);
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::polyBoundaryMeshEntries::polyBoundaryMeshEntries(const IOobject& io)
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
Foam::polyBoundaryMeshEntries::readContents(const IOobject& io)
{
    polyBoundaryMeshEntries reader(io);

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

void Foam::polyBoundaryMeshEntries::removeProcPatches(PtrList<entry>& entries)
{
    // Truncate at the first processor patch entry

    label nNonProcessor = entries.size();

    forAll(entries, patchi)
    {
        const dictionary& dict = entries[patchi].dict();

        const word pType = dict.get<word>("type");
        if (pType == processorPolyPatch::typeName)
        {
            nNonProcessor = patchi;
            break;
        }
    }

    entries.resize(nNonProcessor);
}


bool Foam::polyBoundaryMeshEntries::writeEntries
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


Foam::wordList Foam::polyBoundaryMeshEntries::types
(
    const UPtrList<entry>& entries
)
{
    return extract<word>("type", entries, "patch");
}


Foam::labelList Foam::polyBoundaryMeshEntries::patchStarts
(
    const UPtrList<entry>& entries
)
{
    return extract<label>("startFace", entries, 0);
}


Foam::labelList Foam::polyBoundaryMeshEntries::patchSizes
(
    const UPtrList<entry>& entries
)
{
    return extract<label>("nFaces", entries, 0);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyBoundaryMeshEntries::removeProcPatches()
{
    removeProcPatches(*this);
}


Foam::wordList Foam::polyBoundaryMeshEntries::types() const
{
    return extract<word>("type", *this, "patch");
}


Foam::labelList Foam::polyBoundaryMeshEntries::patchStarts() const
{
    return extract<label>("startFace", *this, 0);
}


Foam::labelList Foam::polyBoundaryMeshEntries::patchSizes() const
{
    return extract<label>("nFaces", *this, 0);
}


void Foam::polyBoundaryMeshEntries::writeEntry(Ostream& os) const
{
    writeEntries(os, *this);
}


void Foam::polyBoundaryMeshEntries::writeEntry
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


bool Foam::polyBoundaryMeshEntries::writeData(Ostream& os) const
{
    return writeEntries(os, *this);
}


bool Foam::polyBoundaryMeshEntries::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    streamOpt.compression(IOstreamOption::UNCOMPRESSED);
    return regIOobject::writeObject(streamOpt, writeOnProc);
}


// ************************************************************************* //
