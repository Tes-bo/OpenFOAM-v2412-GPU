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

#include "LduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::autoPtr<typename Foam::LduMatrix<Type, DType, LUType>::preconditioner>
Foam::LduMatrix<Type, DType, LUType>::preconditioner::New
(
    const solver& sol,
    const dictionary& preconditionerDict
)
{
    const word preconditionerName
    (
        preconditionerDict.get<word>("preconditioner")
    );

    if (sol.matrix().symmetric())
    {
        auto* ctorPtr = symMatrixConstructorTable(preconditionerName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                preconditionerDict,
                "symmetric matrix preconditioner",
                preconditionerName,
                *symMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            ctorPtr
            (
                sol,
                preconditionerDict
            )
        );
    }
    else if (sol.matrix().asymmetric())
    {
        auto* ctorPtr = asymMatrixConstructorTable(preconditionerName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                preconditionerDict,
                "asymmetric matrix preconditioner",
                preconditionerName,
                *asymMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<typename LduMatrix<Type, DType, LUType>::preconditioner>
        (
            ctorPtr
            (
                sol,
                preconditionerDict
            )
        );
    }

    FatalIOErrorInFunction(preconditionerDict)
        << "Cannot precondition incomplete matrix, "
           "no diagonal or off-diagonal coefficient"
        << exit(FatalIOError);

    return nullptr;
}


// ************************************************************************* //
