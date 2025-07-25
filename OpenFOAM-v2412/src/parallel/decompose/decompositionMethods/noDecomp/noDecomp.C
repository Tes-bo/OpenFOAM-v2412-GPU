/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

#include "noDecomp.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeName(noDecomp);
    addNamedToRunTimeSelectionTable
    (
        decompositionMethod,
        noDecomp,
        dictionary,
        none
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::noDecomp::noDecomp(const label numDomains)
:
    decompositionMethod(numDomains)
{}


Foam::noDecomp::noDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::noDecomp::decompose
(
    const pointField& points,
    const scalarField&
) const
{
    return labelList(points.size(), UPstream::myProcNo());
}


Foam::labelList Foam::noDecomp::decompose
(
    const polyMesh& mesh,
    const pointField&,
    const scalarField&
) const
{
    return labelList(mesh.nCells(), UPstream::myProcNo());
}


Foam::labelList Foam::noDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    return labelList(globalCellCells.size(), UPstream::myProcNo());
}


Foam::labelList Foam::noDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    return labelList(globalCellCells.size(), UPstream::myProcNo());
}


// ************************************************************************* //
