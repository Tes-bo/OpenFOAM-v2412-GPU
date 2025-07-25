/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "GAMGPreconditioner.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(GAMGPreconditioner, 0);

    lduMatrix::preconditioner::addsymMatrixConstructorToTable
    <GAMGPreconditioner> addGAMGPreconditionerSymMatrixConstructorToTable_;

    lduMatrix::preconditioner::addasymMatrixConstructorToTable
    <GAMGPreconditioner> addGAMGPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GAMGPreconditioner::GAMGPreconditioner
(
    const lduMatrix::solver& sol,
    const dictionary& solverControls
)
:
    GAMGSolver
    (
        sol.fieldName(),
        sol.matrix(),
        sol.interfaceBouCoeffs(),
        sol.interfaceIntCoeffs(),
        sol.interfaces(),
        solverControls
    ),
    lduMatrix::preconditioner(sol),
    nVcycles_(2)
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GAMGPreconditioner::readControls()
{
    GAMGSolver::readControls();
    nVcycles_ = controlDict_.getOrDefault<label>("nVcycles", 2);
}


void Foam::GAMGPreconditioner::precondition
(
    solveScalarField& wA,
    const solveScalarField& rA_ss,
    const direction cmpt
) const
{
    wA = Zero;
    solveScalarField AwA(wA.size());
    solveScalarField finestCorrection(wA.size());
    solveScalarField finestResidual(rA_ss);

    // Create coarse grid correction fields
    PtrList<solveScalarField> coarseCorrFields;

    // Create coarse grid sources
    PtrList<solveScalarField> coarseSources;

    // Create the smoothers for all levels
    PtrList<lduMatrix::smoother> smoothers;

    // Scratch fields if processor-agglomerated coarse level meshes
    // are bigger than original. Usually not needed
    solveScalarField ApsiScratch;
    solveScalarField finestCorrectionScratch;

    // Initialise the above data structures
    initVcycle
    (
        coarseCorrFields,
        coarseSources,
        smoothers,
        ApsiScratch,
        finestCorrectionScratch
    );

    // Adapt solveScalarField back to scalarField (as required)
    ConstPrecisionAdaptor<scalar, solveScalar> rA_adaptor(rA_ss);
    const scalarField& rA = rA_adaptor.cref();

    for (label cycle=0; cycle<nVcycles_; cycle++)
    {
        Vcycle
        (
            smoothers,
            wA,
            rA,
            AwA,
            finestCorrection,
            finestResidual,

            (ApsiScratch.size() ? ApsiScratch : AwA),
            (
                finestCorrectionScratch.size()
              ? finestCorrectionScratch
              : finestCorrection
            ),

            coarseCorrFields,
            coarseSources,
            cmpt
        );

        if (cycle < nVcycles_-1)
        {
            // Calculate finest level residual field
            matrix_.Amul(AwA, wA, interfaceBouCoeffs_, interfaces_, cmpt);
            finestResidual = rA_ss;
            finestResidual -= AwA;
        }
    }
}


// ************************************************************************* //
