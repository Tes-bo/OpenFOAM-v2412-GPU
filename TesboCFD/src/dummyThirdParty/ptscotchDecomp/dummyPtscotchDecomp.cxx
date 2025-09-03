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

#include "ptscotchDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "Time.H"

static const char* notImplementedMessage =
"Attempted to use <ptscotch> without the ptscotchDecomp library loaded.\n"
"This message is from the dummy ptscotchDecomp stub library instead.\n\n"
"Please install <ptscotch> and ensure libptscotch.so is in LD_LIBRARY_PATH.\n"
"The ptscotchDecomp library can then be built from "
"src/parallel/decompose/ptscotchDecomp.\n"
"Dynamically loading or linking this library will add "
"<ptscotch> as a decomposition method.\n";


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ptscotchDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        ptscotchDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::ptscotchDecomp::decompose
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cWeights,
    labelList& finalDecomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ptscotchDecomp::ptscotchDecomp(const label numDomains)
:
    decompositionMethod(numDomains),
    coeffsDict_()
{}


Foam::ptscotchDecomp::ptscotchDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    decompositionMethod(decompDict, regionName),
    coeffsDict_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::ptscotchDecomp::decompose
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


Foam::labelList Foam::ptscotchDecomp::decompose
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


Foam::labelList Foam::ptscotchDecomp::decompose
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


Foam::labelList Foam::ptscotchDecomp::decompose
(
    const labelListList& globalCellCells,
    const pointField& cellCentres,
    const scalarField& cWeights
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return labelList();
}


// ************************************************************************* //
