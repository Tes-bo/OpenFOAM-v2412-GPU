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

#include "manualRenumber.H"
#include "addToRunTimeSelectionTable.H"
#include "labelIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(manualRenumber);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        manualRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::manualRenumber::manualRenumber(const fileName& file)
:
    renumberMethod(),
    dataFile_(file)
{}


Foam::manualRenumber::manualRenumber(const dictionary& dict)
:
    renumberMethod(dict),
    dataFile_
    (
        dict.optionalSubDict(typeName + "Coeffs").get<fileName>("dataFile")
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::manualRenumber::renumber
(
    const polyMesh& mesh
) const
{
    const label nCells = mesh.nCells();

    labelList newToOld
    (
        labelIOList::readContents
        (
            IOobject
            (
                dataFile_,
                mesh.facesInstance(),
                mesh.thisDb(),
                IOobject::MUST_READ
            )
        )
    );

    // Check if the final renumbering is OK
    if (newToOld.size() != nCells)
    {
        FatalErrorInFunction
            << "Size of renumber list: "
            << newToOld.size() << " != number of cells: " << nCells << nl
            << "Renumbering data read from file " << dataFile_ << endl
            << exit(FatalError);
    }

    // Invert to see if one to one
    labelList oldToNew(nCells, -1);
    forAll(newToOld, i)
    {
        const label origCelli = newToOld[i];

        if (origCelli < 0 || origCelli >= nCells)
        {
            FatalErrorInFunction
                << "Renumbering range error. Index " << i
                << " maps to cell " << origCelli << " from " << nCells << nl
                << "Renumbering data read from file " << dataFile_ << endl
                << exit(FatalError);
        }

        if (oldToNew[origCelli] == -1)
        {
            oldToNew[origCelli] = i;
        }
        else
        {
            FatalErrorInFunction
                << "Renumbering is not one-to-one. Index " << i << " and "
                << oldToNew[origCelli] << " map onto " << origCelli << nl
                << "Renumbering data read from file " << dataFile_ << endl
                << exit(FatalError);
        }
    }

    return newToOld;
}


// ************************************************************************* //
