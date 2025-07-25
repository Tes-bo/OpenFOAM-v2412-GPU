/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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

#include "randomRenumber.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(randomRenumber);

    addToRunTimeSelectionTable
    (
        renumberMethod,
        randomRenumber,
        dictionary
    );
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

labelList randomMap(const label nCells)
{
    Random rndGen(0);

    // Full coverage
    labelList newToOld(Foam::identity(nCells));

    // Fisher-Yates shuffle algorithm
    for (label i = nCells - 1; i > 0; --i)
    {
        label j = rndGen.position<label>(0, i);
        std::swap(newToOld[i], newToOld[j]);
    }

    return newToOld;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::randomRenumber::randomRenumber()
:
    renumberMethod()
{}


Foam::randomRenumber::randomRenumber(const dictionary& dict)
:
    renumberMethod(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::randomRenumber::renumber
(
    const label nCells
) const
{
    return randomMap(nCells);
}


Foam::labelList Foam::randomRenumber::renumber
(
    const pointField& cellCentres
) const
{
    return randomMap(cellCentres.size());
}


Foam::labelList Foam::randomRenumber::renumber
(
    const polyMesh& mesh
) const
{
    return randomMap(mesh.nCells());
}


Foam::labelList Foam::randomRenumber::renumber
(
    const CompactListList<label>& cellCells
) const
{
    return randomMap(cellCells.size());
}


Foam::labelList Foam::randomRenumber::renumber
(
    const labelListList& cellCells
) const
{
    return randomMap(cellCells.size());
}


// ************************************************************************* //
