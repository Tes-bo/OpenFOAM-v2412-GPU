/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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
#include "diagonalSolver.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(lduMatrix::solver, symMatrix);
    defineRunTimeSelectionTable(lduMatrix::solver, asymMatrix);
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::lduMatrix::solver> Foam::lduMatrix::solver::New
(
    const word& solverName,
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
{
    if (matrix.diagonal())
    {
        return autoPtr<lduMatrix::solver>
        (
            new diagonalSolver
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
            )
        );
    }
    else if (matrix.symmetric())
    {
        auto* ctorPtr = symMatrixConstructorTable(solverName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                solverControls,
                "symmetric matrix solver",
                solverName,
                *symMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            ctorPtr
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
            )
        );
    }
    else if (matrix.asymmetric())
    {
        auto* ctorPtr = asymMatrixConstructorTable(solverName);

        if (!ctorPtr)
        {
            FatalIOErrorInLookup
            (
                solverControls,
                "asymmetric matrix solver",
                solverName,
                *asymMatrixConstructorTablePtr_
            ) << exit(FatalIOError);
        }

        return autoPtr<lduMatrix::solver>
        (
            ctorPtr
            (
                fieldName,
                matrix,
                interfaceBouCoeffs,
                interfaceIntCoeffs,
                interfaces,
                solverControls
            )
        );
    }

    FatalIOErrorInFunction(solverControls)
        << "cannot solve incomplete matrix, "
        "no diagonal or off-diagonal coefficient"
        << exit(FatalIOError);

    return nullptr;
}


Foam::autoPtr<Foam::lduMatrix::solver> Foam::lduMatrix::solver::New
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
{
    return New
    (
        solverControls.get<word>("solver"),
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::lduMatrix::solver::solver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    fieldName_(fieldName),
    matrix_(matrix),
    interfaceBouCoeffs_(interfaceBouCoeffs),
    interfaceIntCoeffs_(interfaceIntCoeffs),
    interfaces_(interfaces),
    controlDict_(solverControls),

    log_(1),
    minIter_(0),
    maxIter_(lduMatrix::defaultMaxIter),
    normType_(lduMatrix::normTypes::DEFAULT_NORM),
    tolerance_(lduMatrix::defaultTolerance),
    relTol_(Zero),

    profiling_("lduMatrix::solver." + fieldName)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::lduMatrix::solver::readControls()
{
    log_ = 1;
    minIter_ = 0;
    maxIter_ = lduMatrix::defaultMaxIter;
    normType_ = lduMatrix::normTypes::DEFAULT_NORM;
    tolerance_ = lduMatrix::defaultTolerance;
    relTol_ = 0;

    controlDict_.readIfPresent("log", log_);
    lduMatrix::normTypesNames_.readIfPresent("norm", controlDict_, normType_);
    controlDict_.readIfPresent("minIter", minIter_);
    controlDict_.readIfPresent("maxIter", maxIter_);
    controlDict_.readIfPresent("tolerance", tolerance_);
    controlDict_.readIfPresent("relTol", relTol_);
}


void Foam::lduMatrix::solver::read(const dictionary& solverControls)
{
    controlDict_ = solverControls;
    readControls();
}


Foam::solverPerformance Foam::lduMatrix::solver::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<scalar, solveScalar> tpsi_s(psi);
    return solve
    (
        tpsi_s.ref(),
        ConstPrecisionAdaptor<scalar, solveScalar>(source)(),
        cmpt
    );
}


Foam::solveScalarField::cmptType Foam::lduMatrix::solver::normFactor
(
    const solveScalarField& psi,
    const solveScalarField& source,
    const solveScalarField& Apsi,
    solveScalarField& tmpField,
    const lduMatrix::normTypes normType
) const
{
    switch (normType)
    {
        case lduMatrix::normTypes::NO_NORM :
        {
            break;
        }

        case lduMatrix::normTypes::DEFAULT_NORM :
        case lduMatrix::normTypes::L1_SCALED_NORM :
        {
            // --- Calculate A dot reference value of psi
            matrix_.sumA(tmpField, interfaceBouCoeffs_, interfaces_);

            tmpField *= gAverage(psi, matrix_.mesh().comm());

            return
                gSum
                (
                    (mag(Apsi - tmpField) + mag(source - tmpField))(),
                    matrix_.mesh().comm()
                ) + solverPerformance::small_;

            // Equivalent at convergence:
            // return 2*gSumMag(source) + solverPerformance::small_;
            break;
        }
    }

    // Fall-through: no norm
    return solveScalarField::cmptType(1);
}


// ************************************************************************* //
