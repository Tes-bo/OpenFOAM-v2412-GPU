/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "dlLibraryTable.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TablePtr>
bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry,
    const TablePtr& tablePtr,
    bool verbose
)
{
    List<fileName> libNames;
    dict.readIfPresent(libsEntry, libNames, keyType::LITERAL);

    label nOpen = 0;

    for (const fileName& libName : libNames)
    {
        const label nEntries = (tablePtr ? tablePtr->size() : -1);

        if (dlLibraryTable::open(libName, verbose))
        {
            ++nOpen;

            if (debug && tablePtr != nullptr && tablePtr->size() <= nEntries)
            {
                WarningInFunction
                    << "library " << libName
                    << " did not introduce any new entries"
                    << nl << endl;
            }
        }
    }

    return nOpen && nOpen == libNames.size();
}


// ************************************************************************* //
