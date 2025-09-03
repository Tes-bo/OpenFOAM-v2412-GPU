/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "FDICSmoother.H"
#include "DICPreconditioner.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(FDICSmoother, 0);

    lduMatrix::smoother::addsymMatrixConstructorToTable<FDICSmoother>
        addFDICSmootherSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FDICSmoother::FDICSmoother
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    lduMatrix::smoother
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces
    ),
    rD_(matrix_.diag().size()),
    rDuUpper_(matrix_.upper().size()),
    rDlUpper_(matrix_.upper().size())
{
    solveScalar* __restrict__ rDPtr = rD_.begin();
    solveScalar* __restrict__ rDuUpperPtr = rDuUpper_.begin();
    solveScalar* __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();
    const scalar* const __restrict__ upperPtr =
        matrix_.upper().begin();

    const label nFaces = matrix_.upper().size();

    const scalarField& diag = matrix_.diag();
    std::copy(diag.begin(), diag.end(), rD_.begin());

    DICPreconditioner::calcReciprocalD(rD_, matrix_);

    for (label face=0; face<nFaces; face++)
    {
        rDuUpperPtr[face] = rDPtr[uPtr[face]]*upperPtr[face];
        rDlUpperPtr[face] = rDPtr[lPtr[face]]*upperPtr[face];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FDICSmoother::smooth
(
    solveScalarField& psi,
    const scalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    const solveScalar* const __restrict__ rDuUpperPtr = rDuUpper_.begin();
    const solveScalar* const __restrict__ rDlUpperPtr = rDlUpper_.begin();

    const label* const __restrict__ uPtr =
        matrix_.lduAddr().upperAddr().begin();
    const label* const __restrict__ lPtr =
        matrix_.lduAddr().lowerAddr().begin();

    // Temporary storage for the residual
    solveScalarField rA(rD_.size());
    solveScalar* __restrict__ rAPtr = rA.begin();

    for (label sweep=0; sweep<nSweeps; sweep++)
    {
        matrix_.residual
        (
            rA,
            psi,
            source,
            interfaceBouCoeffs_,
            interfaces_,
            cmpt
        );

        forAll(rA, i)
        {
            rA[i] *= rD_[i];
        }

        const label nFaces = matrix_.upper().size();
        for (label face=0; face<nFaces; face++)
        {
            rAPtr[uPtr[face]] -= rDuUpperPtr[face]*rAPtr[lPtr[face]];
        }

        const label nFacesM1 = nFaces - 1;
        for (label face=nFacesM1; face>=0; face--)
        {
            rAPtr[lPtr[face]] -= rDlUpperPtr[face]*rAPtr[uPtr[face]];
        }

        psi += rA;
    }
}


void Foam::FDICSmoother::scalarSmooth
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const label nSweeps
) const
{
    smooth
    (
        psi,
        ConstPrecisionAdaptor<scalar, solveScalar>(source),
        cmpt,
        nSweeps
    );
}


// ************************************************************************* //
