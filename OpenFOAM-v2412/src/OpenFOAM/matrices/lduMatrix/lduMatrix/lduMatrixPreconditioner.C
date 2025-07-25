/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "lduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduMatrix::preconditioner, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::preconditioner, asymMatrix);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::lduMatrix::preconditioner::getName
(
    const dictionary& solverControls
)
{
    word name;

    // Handle primitive or dictionary entry
    const entry& e =
        solverControls.lookupEntry("preconditioner", keyType::LITERAL);

    if (e.isDict())
    {
        e.dict().readEntry("preconditioner", name);
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::lduMatrix::preconditioner>
Foam::lduMatrix::preconditioner::New
(
    const solver& sol,
    const dictionary& solverControls
)
{
    word name;

    // Handle primitive or dictionary entry

    const entry& e =
        solverControls.lookupEntry("preconditioner", keyType::LITERAL);

    if (e.isDict())
    {
        e.dict().readEntry("preconditioner", name);
    }
    else
    {
        e.stream() >> name;
    }

    const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (sol.matrix().symmetric())
    {
        auto* ctorPtr = symMatrixConstructorTable(name);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                controls,
                "symmetric matrix preconditioner",
                name,
                *symMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            ctorPtr
            (
                sol,
                controls
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        auto* ctorPtr = asymMatrixConstructorTable(name);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                controls,
                "asymmetric matrix preconditioner",
                name,
                *asymMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::preconditioner>
        (
            ctorPtr
            (
                sol,
                controls
            )
        );
    }

    FatalIOErrorInFunction(controls)
        << "cannot solve incomplete matrix, "
           "no diagonal or off-diagonal coefficient"
        << exit(FatalIOError);

    return nullptr;
}


// ************************************************************************* //
