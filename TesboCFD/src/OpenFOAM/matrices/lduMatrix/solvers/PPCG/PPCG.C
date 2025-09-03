/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 Mattijs Janssens
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "PPCG.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PPCG, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<PPCG>
        addPPCGSymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::PPCG::gSumMagProd
(
    FixedList<solveScalar, 3>& globalSum,
    const solveScalarField& a,
    const solveScalarField& b,
    const solveScalarField& c,
    const solveScalarField& sumMag,
    UPstream::Request& request,
    const label comm
)
{
    const label nCells = a.size();

    globalSum = 0.0;
    for (label cell=0; cell<nCells; ++cell)
    {
        globalSum[0] += a[cell]*b[cell];    // sumProd(a, b)
        globalSum[1] += a[cell]*c[cell];    // sumProd(a, c)
        globalSum[2] += mag(sumMag[cell]);
    }

    if (UPstream::parRun())
    {
        Foam::reduce
        (
            globalSum.data(),
            globalSum.size(),
            sumOp<solveScalar>(),
            UPstream::msgType(),  // (ignored): direct MPI call
            comm,
            request
        );
    }
}


Foam::solverPerformance Foam::PPCG::scalarSolveCG
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt,
    const bool cgMode
) const
{
    // --- Setup class containing solver performance data
    solverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(controlDict_) + type(),
        fieldName_
    );

    const label comm = matrix().mesh().comm();
    const label nCells = psi.size();
    solveScalarField w(nCells);

    // --- Calculate A.psi
    matrix_.Amul(w, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    // --- Calculate initial residual field
    solveScalarField r(source - w);

    // --- Calculate normalisation factor
    solveScalarField p(nCells);
    const solveScalar normFactor = this->normFactor(psi, source, w, p);

    if ((log_ >= 2) || (lduMatrix::debug >= 2))
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // --- Select and construct the preconditioner
    if (!preconPtr_)
    {
        preconPtr_ = lduMatrix::preconditioner::New
        (
            *this,
            controlDict_
        );
    }

    // --- Precondition residual (= u0)
    solveScalarField u(nCells);
    preconPtr_->precondition(u, r, cmpt);

    // --- Calculate A*u - reuse w
    matrix_.Amul(w, u, interfaceBouCoeffs_, interfaces_, cmpt);


    // State
    solveScalarField s(nCells);
    solveScalarField q(nCells);
    solveScalarField z(nCells);

    solveScalarField m(nCells);

    FixedList<solveScalar, 3> globalSum;
    UPstream::Request outstandingRequest;
    if (cgMode)
    {
        // --- Start global reductions for inner products
        gSumMagProd(globalSum, u, r, w, r, outstandingRequest, comm);

        // --- Precondition residual
        preconPtr_->precondition(m, w, cmpt);
    }
    else
    {
        // --- Precondition residual
        preconPtr_->precondition(m, w, cmpt);

        // --- Start global reductions for inner products
        gSumMagProd(globalSum, w, u, m, r, outstandingRequest, comm);
    }

    // --- Calculate A*m
    solveScalarField n(nCells);
    matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);

    solveScalar alpha = 0.0;
    solveScalar gamma = 0.0;

    // --- Solver iteration
    for
    (
        solverPerf.nIterations() = 0;
        solverPerf.nIterations() < maxIter_;
        solverPerf.nIterations()++
    )
    {
        // Make sure gamma,delta are available
        outstandingRequest.wait();

        const solveScalar gammaOld = gamma;
        gamma = globalSum[0];
        const solveScalar delta = globalSum[1];

        solverPerf.finalResidual() = globalSum[2]/normFactor;
        if (solverPerf.nIterations() == 0)
        {
            solverPerf.initialResidual() = solverPerf.finalResidual();
        }

        // Check convergence (bypass if not enough iterations yet)
        if
        (
            (minIter_ <= 0 || solverPerf.nIterations() >= minIter_)
         && solverPerf.checkConvergence(tolerance_, relTol_, log_)
        )
        {
            break;
        }


        if (solverPerf.nIterations() == 0)
        {
            alpha = gamma/delta;
            z = n;
            q = m;
            s = w;
            p = u;
        }
        else
        {
            const solveScalar beta = gamma/gammaOld;
            alpha = gamma/(delta-beta*gamma/alpha);

            for (label cell=0; cell<nCells; ++cell)
            {
                z[cell] = n[cell] + beta*z[cell];
                q[cell] = m[cell] + beta*q[cell];
                s[cell] = w[cell] + beta*s[cell];
                p[cell] = u[cell] + beta*p[cell];
            }
        }

        for (label cell=0; cell<nCells; ++cell)
        {
            psi[cell] += alpha*p[cell];
            r[cell] -= alpha*s[cell];
            u[cell] -= alpha*q[cell];
            w[cell] -= alpha*z[cell];
        }

        if (cgMode)
        {
            // --- Start global reductions for inner products
            gSumMagProd(globalSum, u, r, w, r, outstandingRequest, comm);

            // --- Precondition residual
            preconPtr_->precondition(m, w, cmpt);
        }
        else
        {
            // --- Precondition residual
            preconPtr_->precondition(m, w, cmpt);

            // --- Start global reductions for inner products
            gSumMagProd(globalSum, w, u, m, r, outstandingRequest, comm);
        }

        // --- Calculate A*m
        matrix_.Amul(n, m, interfaceBouCoeffs_, interfaces_, cmpt);
    }

    // Cleanup any outstanding requests
    outstandingRequest.wait();

    if (preconPtr_)
    {
        preconPtr_->setFinished(solverPerf);
    }

    //TBD
    //matrix().setResidualField
    //(
    //    ConstPrecisionAdaptor<scalar, solveScalar>(rA)(),
    //    fieldName_,
    //    false
    //);

    return solverPerf;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PPCG::PPCG
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::PPCG::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolveCG
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt,
        true   // operate in conjugate-gradient mode
    );
}


// ************************************************************************* //
