/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "decomposedBlockData.H"
#include "dictionary.H"
#include "foamVersion.H"
#include "objectRegistry.H"
#include "SpanStream.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Like Ostream::writeEntry, but with fewer spaces
template<class T>
static inline void writeHeaderEntry
(
    Ostream& os,
    const word& key,
    const T& value
)
{
    os.indent();
    os.write(key);

    label padding = (12 - label(key.size()));

    // Write padding spaces (always at least one)
    do
    {
        os.write(char(token::SPACE));
    }
    while (--padding > 0);

    os << value << char(token::END_STATEMENT) << nl;
}

} // End namespace Foam


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::decomposedBlockData::writeHeaderContent
(
    Ostream& os,
    IOstreamOption streamOptContainer,
    const word& objectType,
    const string& note,
    const fileName& location,
    const word& objectName
)
{
    // Standard header entries
    writeHeaderEntry(os, "version", streamOptContainer.version());
    writeHeaderEntry(os, "format", streamOptContainer.format());
    writeHeaderEntry(os, "arch", foamVersion::buildArch);

    if (!note.empty())
    {
        writeHeaderEntry(os, "note", note);
    }

    if (objectType.empty())
    {
        // Empty type not allowed - use 'dictionary' fallback
        writeHeaderEntry(os, "class", word("dictionary"));
    }
    else
    {
        writeHeaderEntry(os, "class", objectType);
    }

    if (!location.empty())
    {
        writeHeaderEntry(os, "location", location);
    }
    writeHeaderEntry(os, "object", objectName);
}


// * * * * * * * * * * * * * * * Members Functions * * * * * * * * * * * * * //

bool Foam::decomposedBlockData::readHeader(IOobject& io, Istream& is)
{
    dictionary headerDict;

    // Read the regular "FoamFile" header
    bool ok = io.readHeader(headerDict, is);

    if (decomposedBlockData::isCollatedType(io))
    {
        // Quick information - extract from "data.class"
        if (headerDict.readIfPresent("data.class", io.headerClassName()))
        {
            return ok;
        }

        {
            // Master-only reading of header
            List<char> charData;
            decomposedBlockData::readBlockEntry(is, charData);

            ISpanStream headerStream(charData);
            headerStream.name() = is.name();

            ok = io.readHeader(headerStream);
        }
    }

    return ok;
}


void Foam::decomposedBlockData::writeHeader
(
    Ostream& os,
    IOstreamOption streamOptContainer,
    const word& objectType,
    const string& note,
    const fileName& location,
    const word& objectName,
    const dictionary& extraEntries
)
{
    if (IOobject::bannerEnabled())
    {
        IOobject::writeBanner(os);
    }

    os.beginBlock("FoamFile");

    decomposedBlockData::writeHeaderContent
    (
        os,
        streamOptContainer,
        objectType,
        note,
        location,
        objectName
    );

    if (!extraEntries.empty())
    {
        extraEntries.writeEntries(os);
    }

    os.endBlock();

    if (IOobject::bannerEnabled())
    {
        IOobject::writeDivider(os) << nl;
    }
}


void Foam::decomposedBlockData::writeExtraHeaderContent
(
    dictionary& dict,
    IOstreamOption streamOptData,
    const IOobject& io
)
{
    dict.set("data.format", streamOptData.format());
    dict.set("data.class", io.type());

    // Deep-copy of meta-data (if any)
    const dictionary* metaDataDict = io.findMetaData();
    if (metaDataDict && !metaDataDict->empty())
    {
        dict.add("meta", *metaDataDict);
    }
}


void Foam::decomposedBlockData::writeHeader
(
    Ostream& os,
    IOstreamOption streamOptData,
    const IOobject& io
)
{
    if (IOobject::bannerEnabled())
    {
        IOobject::writeBanner(os);
    }

    os.beginBlock("FoamFile");

    decomposedBlockData::writeHeaderContent
    (
        os,
        static_cast<IOstreamOption>(os),    // streamOpt container
        decomposedBlockData::typeName,      // class
        io.note(),
        (io.instance()/io.db().dbDir()/io.local()),  // location
        io.name()
    );

    {
        writeHeaderEntry(os, "data.format", streamOptData.format());
        writeHeaderEntry(os, "data.class", io.type());
    }

    // Meta-data (if any)
    const dictionary* metaDataDict = io.findMetaData();
    if (metaDataDict && !metaDataDict->empty())
    {
        metaDataDict->writeEntry("meta", os);
    }

    os.endBlock();

    if (IOobject::bannerEnabled())
    {
        IOobject::writeDivider(os) << nl;
    }
}


// ************************************************************************* //
