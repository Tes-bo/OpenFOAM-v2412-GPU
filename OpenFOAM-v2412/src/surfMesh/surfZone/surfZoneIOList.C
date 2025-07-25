/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(surfZoneIOList);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::surfZoneIOList::readIOcontents()
{
    if
    (
        this->isReadRequired()
    )
    {
        surfZoneList& zones = *this;

        // Read entries
        Istream& is = readStream(typeName);

        PtrList<entry> entries(is);
        zones.resize(entries.size());

        // Transcribe
        label startOffset = 0;
        forAll(zones, zonei)
        {
            zones[zonei] = surfZone
            (
                entries[zonei].keyword(),
                entries[zonei].dict(),
                zonei
            );

            if (zones[zonei].start() != startOffset)
            {
                FatalErrorInFunction
                    << "surfZones are not ordered. Start of zone " << zonei
                    << " does not correspond to sum of preceding zones." << nl
                    << "while reading " << this->objectRelPath() << endl
                    << exit(FatalError);
            }

            startOffset += zones[zonei].size();
        }

        is.check(FUNCTION_NAME);
        close();
        return true;
    }

    // Nothing read
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io
)
:
    regIOobject(io),
    surfZoneList()
{
    readIOcontents();  // allowOptionalRead = false
}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    const UList<surfZone>& content
)
:
    regIOobject(io),
    surfZoneList(content)
{}


Foam::surfZoneIOList::surfZoneIOList
(
    const IOobject& io,
    surfZoneList&& content
)
:
    regIOobject(io),
    surfZoneList(std::move(content))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::surfZoneIOList::writeData(Ostream& os) const
{
    const surfZoneList& zones = *this;
    const label len = zones.size();

    if (len)
    {
        os  << len << nl << token::BEGIN_LIST << incrIndent << nl;

        for (const surfZone& zn : zones)
        {
            zn.write(os);
        }

        os  << decrIndent << token::END_LIST;
    }
    else
    {
        os  << len << token::BEGIN_LIST << token::END_LIST;
    }

    return os.good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::surfZoneIOList::operator=(const surfZoneIOList& rhs)
{
    surfZoneList::operator=(rhs);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIOList& zones)
{
    zones.writeData(os);

    return os;
}


// ************************************************************************* //
