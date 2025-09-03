/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "scotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"Attempted to use <scotch> without the scotchDecomp library loaded.\n"
"This message is from the dummy scotchDecomp stub library instead.\n\n"
"Please install <scotch> and ensure libscotch.so is in LD_LIBRARY_PATH.\n"
"The scotchDecomp library can then be built from "
"src/parallel/decompose/scotchDecomp.\n"
"Dynamically loading or linking this library will add "
"<scotch> as a decomposition method.\n";


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(scotchDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        scotchDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::scotchDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& decomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << nl
        << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scotchDecomp::scotchDecomp(const label numDomains)
:
    metisLikeDecomp(numDomains)
{}


Foam::scotchDecomp::scotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp("scotch", decompDict, regionName)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const pointField& points,
    const scalarField& pointWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const polyMesh& mesh,
    const labelList& agglom,
    const pointField& agglomPoints,
    const scalarField& agglomWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const CompactListList<label>& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList();
}


Foam::labelList Foam::scotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField&,
    const scalarField&
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList();
}


// ************************************************************************* //
