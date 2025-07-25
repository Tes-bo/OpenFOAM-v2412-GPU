/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "metisDecomp.H"
#include "addToRunTimeSelectionTable.H"

static const char* notImplementedMessage =
"Attempted to use <metis> without the metisDecomp library loaded.\n"
"This message is from the dummy metisDecomp stub library instead.\n\n"
"Please install <metis> and ensure libmetis.so is in LD_LIBRARY_PATH.\n"
"The metisDecomp library can then be built from "
"src/parallel/decompose/metisDecomp.\n"
"Dynamically loading or linking this library will add "
"<metis> as a decomposition method.\n";

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(metisDecomp, 0);
    addToRunTimeSelectionTable
    (
        decompositionMethod,
        metisDecomp,
        dictionary
    );
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

Foam::label Foam::metisDecomp::decomposeSerial
(
    const labelList& adjncy,
    const labelList& xadj,
    const List<scalar>& cellWeights,
    labelList& decomp
) const
{
    FatalErrorInFunction
        << notImplementedMessage << exit(FatalError);

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::metisDecomp::metisDecomp
(
    const dictionary& decompDict,
    const word& regionName
)
:
    metisLikeDecomp("metis", decompDict, regionName)
{}


// ************************************************************************* //
