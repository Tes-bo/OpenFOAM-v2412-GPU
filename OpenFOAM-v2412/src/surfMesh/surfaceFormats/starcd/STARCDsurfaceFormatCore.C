/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "STARCDsurfaceFormatCore.H"
#include "clock.H"
#include "regExp.H"
#include "IFstream.H"
#include "SubStrings.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// parse things like this:
//     CTNAME  1  someName
// don't bother with the older comma-delimited format

Foam::Map<Foam::word>
Foam::fileFormats::STARCDsurfaceFormatCore::readInpCellTable(ISstream& is)
{
    Map<word> lookup;

    if (!is.good())
    {
        return lookup;
    }

    const regExp ctname
    (
        " *CTNA[^ ]*"        // keyword - min 4 chars
        "[[:space:]]+"       // space delimited
        "([0-9]+)"           // 1: <digits>
        "[[:space:]]+"       // space delimited
        "([^,[:space:]].*)", // 2: <name>
        true                 // ignore case
    );

    string line;
    regExp::results_type groups;

    while (is.good() && is.getLine(line).good())
    {
        if (ctname.match(line, groups))
        {
            const label tableId = readLabel(groups.str(1));
            const word tableName = word::validate(groups.str(2), true);

            if (!tableName.empty())
            {
                lookup.insert(tableId, tableName);
            }
        }
    }

    return lookup;
}


void Foam::fileFormats::STARCDsurfaceFormatCore::writeCase
(
    Ostream& os,
    const UList<point>& pts,
    const label nFaces,
    const UList<surfZone>& zoneLst
)
{
    const word caseName = os.name().stem();

    os  << "! STARCD file written " << clock::dateTime().c_str() << nl
        << "! " << pts.size() << " points, " << nFaces << " faces" << nl
        << "! case " << caseName << nl
        << "! ------------------------------" << nl;

    forAll(zoneLst, zoneI)
    {
        os  << "ctable " << zoneI + 1 << " shell" << " ,,,,,," << nl
            << "ctname " << zoneI + 1 << " "
            << zoneLst[zoneI].name() << nl;
    }

    os  << "! ------------------------------" << nl
        << "*set icvo mxv - 1" << nl
        << "vread " << caseName << ".vrt icvo,,,coded" << nl
        << "cread " << caseName << ".cel icvo,,,add,coded" << nl
        << "*set icvo" << nl
        << "! end" << nl;

    os.flush();
}


// ************************************************************************* //
