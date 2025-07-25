/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "topoSetPointSource.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(topoSetPointSource, word);
    defineRunTimeSelectionTable(topoSetPointSource, istream);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::topoSetPointSource::topoSetPointSource(const polyMesh& mesh)
:
    topoSetSource(mesh)
{}


Foam::topoSetPointSource::topoSetPointSource
(
    const polyMesh& mesh,
    const dictionary& dict
)
:
    topoSetSource(mesh, dict)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::topoSetPointSource>
Foam::topoSetPointSource::New
(
    const word& sourceType,
    const polyMesh& mesh,
    const dictionary& dict
)
{
    auto* ctorPtr = wordConstructorTable(sourceType);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            dict,
            "pointSetSource",
            sourceType,
            *wordConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return autoPtr<topoSetPointSource>(ctorPtr(mesh, dict));
}


Foam::autoPtr<Foam::topoSetPointSource> Foam::topoSetPointSource::New
(
    const word& sourceType,
    const polyMesh& mesh,
    Istream& is
)
{
    auto* ctorPtr = istreamConstructorTable(sourceType);

    if (!ctorPtr)
    {
        FatalErrorInLookup
        (
            "pointSetSource",
            sourceType,
            *istreamConstructorTablePtr_
        ) << exit(FatalError);
    }

    return autoPtr<topoSetPointSource>(ctorPtr(mesh, is));
}


// ************************************************************************* //
