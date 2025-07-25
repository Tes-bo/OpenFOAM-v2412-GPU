/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2013 OpenFOAM Foundation
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

#include "wordIOList.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineCompoundTypeName(List<word>, wordList);
    addCompoundToRunTimeSelectionTable(List<word>, wordList);

    defineTemplateTypeNameAndDebugWithName(wordIOList, "wordList", 0);
    defineTemplateTypeNameAndDebugWithName(wordListIOList, "wordListList", 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::printTable
(
    const UList<wordList>& tbl,
    List<std::string::size_type>& columnWidths,
    Ostream& os,
    bool headerSeparator
)
{
    if (tbl.empty())
    {
        return os;
    }

    // Find maximum width for each column
    columnWidths.resize(tbl.first().size(), std::string::size_type(0));

    forAll(columnWidths, coli)
    {
        auto& colWidth = columnWidths[coli];

        for (const wordList& tblRow : tbl)
        {
            colWidth =
                std::max
                (
                    colWidth,
                    string::size_type(tblRow[coli].length())
                );
        }
    }

    // Print the rows adding spacing for the columns
    for (const wordList& tblRow : tbl)
    {
        forAll(tblRow, coli)
        {
            os  << tblRow[coli];
            for
            (
                string::size_type space = 0;
                space < columnWidths[coli] - tblRow[coli].length() + 2;
                ++space
            )
            {
                os  << ' ';
            }
        }
        os  << nl;

        if (headerSeparator) os  << nl;
        headerSeparator = false;
    }

    return os;
}


Foam::Ostream& Foam::printTable
(
    const UList<wordList>& tbl,
    Ostream& os,
    bool headerSeparator
)
{
    List<std::string::size_type> columnWidths;
    printTable(tbl, columnWidths, os, headerSeparator);
    return os;
}


// ************************************************************************* //
