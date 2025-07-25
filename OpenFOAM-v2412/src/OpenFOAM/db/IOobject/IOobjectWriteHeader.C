/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "IOobject.H"
#include "dictionary.H"
#include "objectRegistry.H"
#include "foamVersion.H"

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


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

// A banner corresponding to this:
//
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  VERSION                               |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

Foam::Ostream&
Foam::IOobject::writeBanner(Ostream& os, const bool noSyntaxHint)
{
    // The version padded with spaces to fit after "Version:  "
    // - initialized with zero-length string to detect if it has been populated
    static char paddedVersion[39] = "";

    if (!*paddedVersion)
    {
        // Populate: like strncpy but without trailing '\0'

        const std::string apiValue(std::to_string(Foam::foamVersion::api));

        std::size_t len = apiValue.length();
        if (len > 38)
        {
            len = 38;
        }

        std::memset(paddedVersion, ' ', 38);
        std::memcpy(paddedVersion, apiValue.c_str(), len);
        paddedVersion[38] = '\0';
    }

    os  <<
        "/*--------------------------------";

    if (noSyntaxHint)
    {
        // Without syntax hint
        os  << "---------";
    }
    else
    {
        // With syntax hint
        os  << "*- C++ -*";
    }

    os  <<
        "----------------------------------*\\\n"
        "| =========                 |"
        "                                                 |\n"
        "| \\\\      /  F ield         |"
        " OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     |"
        " Version:  " << paddedVersion << "|\n"
        "|   \\\\  /    A nd           |"
        " Website:  www.openfoam.com                      |\n"
        "|    \\\\/     M anipulation  |"
        "                                                 |\n"
        "\\*-----------------------------------------"
        "----------------------------------*/\n";

    return os;
}


Foam::Ostream& Foam::IOobject::writeDivider(Ostream& os)
{
    os  <<
        "// * * * * * * * * * * * * * * * * * "
        "* * * * * * * * * * * * * * * * * * * * //\n";

    return os;
}


Foam::Ostream& Foam::IOobject::writeEndDivider(Ostream& os)
{
    os  << "\n\n"
        "// *****************************************"
        "******************************** //\n";

    return os;
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::IOobject::writeHeaderContent
(
    Ostream& os,
    const IOobject& io,
    const word& objectType,
    const dictionary* metaDataDict
)
{
    // Standard header entries
    writeHeaderEntry(os, "version", os.version());
    writeHeaderEntry(os, "format", os.format());
    writeHeaderEntry(os, "arch", foamVersion::buildArch);

    if (!io.note().empty())
    {
        writeHeaderEntry(os, "note", io.note());
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

    writeHeaderEntry(os, "location", io.instance()/io.db().dbDir()/io.local());
    writeHeaderEntry(os, "object", io.name());

    // Meta-data (if any)
    if (metaDataDict && !metaDataDict->empty())
    {
        metaDataDict->writeEntry("meta", os);
    }
}


void Foam::IOobject::writeHeaderContent
(
    dictionary& dict,
    const IOobject& io,
    const word& objectType,
    IOstreamOption streamOpt,
    const dictionary* metaDataDict
)
{
    // Standard header entries
    dict.set("version", streamOpt.version());
    dict.set("format", streamOpt.format());
    dict.set("arch", foamVersion::buildArch);

    if (!io.note().empty())
    {
        dict.set("note", io.note());
    }

    if (objectType.empty())
    {
        // Empty type not allowed - use 'dictionary' fallback
        dict.set("class", word("dictionary"));
    }
    else
    {
        dict.set("class", objectType);
    }

    dict.set("location", io.instance()/io.db().dbDir()/io.local());
    dict.set("object", io.name());

    // Deep-copy of meta-data (if any)
    if (metaDataDict && !metaDataDict->empty())
    {
        dict.add("meta", *metaDataDict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOobject::writeHeader
(
    Ostream& os,
    const word& objectType
) const
{
    if (!os.good())
    {
        InfoInFunction
            << "No stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    if (IOobject::bannerEnabled())
    {
        IOobject::writeBanner(os);
    }

    os.beginBlock("FoamFile");

    // Standard header entries
    IOobject::writeHeaderContent
    (
        os,
        *this,
        objectType,
        this->findMetaData()
    );

    os.endBlock();

    if (IOobject::bannerEnabled())
    {
        IOobject::writeDivider(os) << nl;
    }

    return true;
}


bool Foam::IOobject::writeHeader(Ostream& os) const
{
    return IOobject::writeHeader(os, this->type());
}


void Foam::IOobject::writeHeader
(
    dictionary& dict,
    const word& objectType,
    IOstreamOption streamOpt
) const
{
    IOobject::writeHeaderContent
    (
        dict,
        *this,
        objectType,
        streamOpt,
        this->findMetaData()
    );
}


void Foam::IOobject::writeHeader
(
    dictionary& dict,
    IOstreamOption streamOpt
) const
{
    IOobject::writeHeader(dict, this->type(), streamOpt);
}


// ************************************************************************* //
