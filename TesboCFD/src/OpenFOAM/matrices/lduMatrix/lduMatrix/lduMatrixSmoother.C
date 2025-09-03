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
    defineRunTimeSelectionTable(lduMatrix::smoother, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::smoother, asymMatrix);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word
Foam::lduMatrix::smoother::getName
(
    const dictionary& solverControls
)
{
    word name;

    // Handle primitive or dictionary entry
    const entry& e =
        solverControls.lookupEntry("smoother", keyType::LITERAL);

    if (e.isDict())
    {
        e.dict().readEntry("smoother", name);
    }
    else
    {
        e.stream() >> name;
    }

    return name;
}


Foam::autoPtr<Foam::lduMatrix::smoother> Foam::lduMatrix::smoother::New
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
{
    word name;

    // Handle primitive or dictionary entry
    const entry& e =
        solverControls.lookupEntry("smoother", keyType::LITERAL);

    if (e.isDict())
    {
        e.dict().readEntry("smoother", name);
    }
    else
    {
        e.stream() >> name;
    }

    // not (yet?) needed:
    // const dictionary& controls = e.isDict() ? e.dict() : dictionary::null;

    if (matrix.symmetric())
    {
        auto* ctorPtr = symMatrixConstructorTable(name);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                solverControls,
                "symmetric matrix smoother",
                name,
                *symMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            ctorPtr
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }
    else if (matrix.asymmetric())
    {
        auto* ctorPtr = asymMatrixConstructorTable(name);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                solverControls,
                "asymmetric matrix smoother",
                name,
                *asymMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::smoother>
        (
            ctorPtr
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces
            )
        );
    }

    FatalIOErrorInFunction(solverControls)
        << "cannot solve incomplete matrix, "
        "no diagonal or off-diagonal coefficient"
        << exit(FatalIOError);

    return nullptr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::smoother::smoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    fieldName_(fieldName),
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaceIntCoeffs_(interfaceIntCoeffs),
    interfaces_(interfaces)
{}


// ************************************************************************* //
