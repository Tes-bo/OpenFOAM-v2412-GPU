/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "scalarRanges.H"
#include "stringOps.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::scalarRanges Foam::scalarRanges::parse
(
    const std::string& str,
    bool report
)
{
    const auto items = stringOps::splitAny(str, " ,;");

    scalarRanges ranges(items.size());

    label n = 0;

    for (const auto& item : items)
    {
        const std::string s(item.str());

        scalarRange& range = ranges[n];

        if (scalarRange::parse(s, range))
        {
            ++n;
        }
        else if (report)
        {
            InfoErr<< "Bad scalar-range parsing: " << s << endl;
        }
    }

    ranges.resize(n);

    return ranges;
}


// ************************************************************************* //
