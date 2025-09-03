/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "noRenumber.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(noRenumber);

    addNamedToRunTimeSelectionTable
    (
        renumberMethod,
        noRenumber,
        dictionary,
        none
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noRenumber::noRenumber()
:
    renumberMethod()
{}


Foam::noRenumber::noRenumber(const dictionary& dict)
:
    renumberMethod(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::noRenumber::renumber
(
    const label nCells
) const
{
    return Foam::identity(nCells);
}


Foam::labelList Foam::noRenumber::renumber
(
    const pointField& cellCentres
) const
{
    return Foam::identity(cellCentres.size());
}


Foam::labelList Foam::noRenumber::renumber
(
    const polyMesh& mesh
) const
{
    return Foam::identity(mesh.nCells());
}


Foam::labelList Foam::noRenumber::renumber
(
    const CompactListList<label>& cellCells
) const
{
    return Foam::identity(cellCells.size());
}


Foam::labelList Foam::noRenumber::renumber
(
    const labelListList& cellCells
) const
{
    return Foam::identity(cellCells.size());
}


// ************************************************************************* //
