/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "CuthillMcKeeRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "bandCompression.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(CuthillMcKeeRenumber);
    defineTypeName(reverseCuthillMcKeeRenumber);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        CuthillMcKeeRenumber,
        dictionary
    );

    addToRunTimeSelectionTable
    (
        renumberMethod,
        reverseCuthillMcKeeRenumber,
        dictionary
    );

    // Select under the name "RCM"
    addNamedToRunTimeSelectionTable
    (
        renumberMethod,
        reverseCuthillMcKeeRenumber,
        dictionary,
        RCM
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CuthillMcKeeRenumber::CuthillMcKeeRenumber(const bool reverse)
:
    renumberMethod(),
    reverse_(reverse)
{}


Foam::CuthillMcKeeRenumber::CuthillMcKeeRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    reverse_
    (
        dict.optionalSubDict(typeName + "Coeffs")
            .getOrDefault("reverse", false)
    )
{}


Foam::CuthillMcKeeRenumber::CuthillMcKeeRenumber
(
    const dictionary& dict,
    const bool reverse
)
:
    renumberMethod(dict),
    reverse_(reverse)
{}


Foam::reverseCuthillMcKeeRenumber::reverseCuthillMcKeeRenumber()
:
    CuthillMcKeeRenumber(true)  // reverse = true
{}


Foam::reverseCuthillMcKeeRenumber::reverseCuthillMcKeeRenumber
(
    const dictionary& dict
)
:
    CuthillMcKeeRenumber(dict, true)  // reverse = true
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const polyMesh& mesh
) const
{
    labelList orderedToOld = meshTools::bandCompression(mesh);

    if (reverse_)
    {
        Foam::reverse(orderedToOld);
    }

    return orderedToOld;
}


Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const CompactListList<label>& cellCells
) const
{
    labelList orderedToOld = meshTools::bandCompression(cellCells);

    if (reverse_)
    {
        Foam::reverse(orderedToOld);
    }

    return orderedToOld;
}


Foam::labelList Foam::CuthillMcKeeRenumber::renumber
(
    const labelListList& cellCells
) const
{
    labelList orderedToOld = meshTools::bandCompression(cellCells);

    if (reverse_)
    {
        Foam::reverse(orderedToOld);
    }

    return orderedToOld;
}


// Foam::labelList Foam::CuthillMcKeeRenumber::renumber
// (
//     const labelUList& cellCells,
//     const labelUList& offsets
// ) const
// {
//     labelList orderedToOld = meshTools::bandCompression(cellCells, offsets);
//
//     if (reverse_)
//     {
//         Foam::reverse(orderedToOld);
//     }
//
//     return orderedToOld;
// }


// ************************************************************************* //
