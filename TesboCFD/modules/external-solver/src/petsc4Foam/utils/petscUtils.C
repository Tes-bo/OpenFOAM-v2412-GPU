/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
    Copyright (C) 2019 Simone Bna
    Copyright (C) 2020 Stefano Zampini
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
#include "error.H"
#include "petscUtils.H"
#include "petscErrorHandling.H"
#include "petscLinearSolverContext.H"

// For older PETSc. Not strictly correct for 64-bit compilations,
// but adequate for transitional code
#ifndef PetscInt_FMT
#define PetscInt_FMT "D"
#endif

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

Foam::solveScalar Foam::gAverage(Vec input)
{
    PetscInt len;
    AssertPETSc(VecGetSize(input, &len));

    if (len)
    {
        PetscScalar val;
        AssertPETSc(VecSum(input, &val));

        return val/len;
    }

    WarningInFunction
        << "Empty PETSc Vec, returning zero" << endl;

    return 0;
}


Foam::solveScalar Foam::gSum(Vec input)
{
    PetscScalar val;
    AssertPETSc(VecSum(input, &val));

    return val;
}


Foam::scalar Foam::PetscUtils::normFactor
(
    Vec AdotPsi,
    Vec psi,
    Vec source,
    Vec ArowsSum
)
{
    // Equivalent to the OpenFOAM normFactor function
    //
    // stabilise
    // (
    //   gSum(cmptMag(Apsi - tmpField) + cmptMag(matrix_.source() - tmpField)),
    //   SolverPerformance<Type>::small_
    // )

    PetscScalar avgPsi;
    {
        AssertPETSc(VecSum(psi, &avgPsi));

        PetscInt len;
        AssertPETSc(VecGetSize(psi, &len));
        avgPsi /= len;
    }

    // TODO: Use Vec primitives?
    const PetscScalar* ArowsSumVecValues;
    AssertPETSc(VecGetArrayRead(ArowsSum, &ArowsSumVecValues));

    const PetscScalar* AdotPsiValues;
    AssertPETSc(VecGetArrayRead(AdotPsi, &AdotPsiValues));

    const PetscScalar* sourceValues;
    AssertPETSc(VecGetArrayRead(source, &sourceValues));

    scalar normFactor{0};

    PetscInt len;
    AssertPETSc(VecGetLocalSize(psi, &len));

    for (PetscInt i=0; i < len; ++i)
    {
        const PetscScalar psiRow = (ArowsSumVecValues[i] * avgPsi);

        normFactor +=
        (
            Foam::mag(AdotPsiValues[i] - psiRow)
          + Foam::mag(sourceValues[i] - psiRow)
        );
    }

    // Restore
    AssertPETSc(VecRestoreArrayRead(ArowsSum, &ArowsSumVecValues));
    AssertPETSc(VecRestoreArrayRead(AdotPsi, &AdotPsiValues));
    AssertPETSc(VecRestoreArrayRead(source, &sourceValues));

    return stabilise
    (
        returnReduce(normFactor, sumOp<scalar>()),
        SolverPerformance<scalar>::small_
    );
}

PetscErrorCode Foam::PetscUtils::foamKSPMonitorFoam
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    void *cctx
)
{
    PetscViewer viewer;
    PetscInt tablevel;
    const char *prefix;
    PetscReal fnorm;
    KSPNormType ntype;

    PetscFunctionBeginUser;
    auto* ctx = static_cast<petscLinearSolverContext*>(cctx);

    // compute L1 norm and rescale by normFactor
    PetscCall(KSPBuildResidual(ksp, ctx->res_l1_w[0], ctx->res_l1_w[1], &ctx->res_l1_w[1]));
    PetscCall(VecNorm(ctx->res_l1_w[1], NORM_1, &fnorm));
    fnorm /= ctx->normFactor;
    PetscCall(PetscViewerASCIIGetStdout(PetscObjectComm((PetscObject)ksp), &viewer));
    PetscCall(PetscObjectGetTabLevel((PetscObject)ksp, &tablevel));
    PetscCall(PetscViewerASCIIAddTab(viewer, tablevel));
    PetscCall(KSPGetOptionsPrefix(ksp, &prefix));
    if (it == 0)
    {
       PetscCall(PetscViewerASCIIPrintf(viewer, "  Residual norms for %s solve.\n", prefix));
    }
    PetscCall(KSPGetNormType(ksp, &ntype));
    if (ntype != KSP_NORM_NONE) // Print both norms, KSP built-in and OpenFOAM built-in
    {
       char normtype[256];

       PetscCall(PetscStrncpy(normtype,KSPNormTypes[ntype],sizeof(normtype)));
       PetscCall(PetscStrtolower(normtype));
       if (ctx->useFoamTest) // we are using foam convergence testing, list foam norm first
       {
           PetscCall(PetscViewerASCIIPrintf
           (
               viewer,
               "%4" PetscInt_FMT " KSP Residual foam norm %14.12e (PETSc %s norm %14.12e)\n",
               it,
               static_cast<double>(fnorm),
               normtype,
               static_cast<double>(rnorm)
           ));
       }
       else // we are using KSP default convergence testing, list KSP norm first
       {
           PetscCall(PetscViewerASCIIPrintf
           (
               viewer,
               "%4" PetscInt_FMT " KSP Residual %s norm %14.12e (foam norm %14.12e)\n",
               it,
               normtype,
               static_cast<double>(rnorm),
               static_cast<double>(fnorm)
           ));
       }
    }
    else // KSP has no norm, list foam norm only
    {
        PetscCall(PetscViewerASCIIPrintf
        (
            viewer,
            "%4" PetscInt_FMT " KSP Residual foam norm %14.12e\n",
            it,
            static_cast<double>(fnorm)
        ));
    }
    PetscCall(PetscViewerASCIISubtractTab(viewer, tablevel));
    PetscFunctionReturn(0);
}


PetscErrorCode Foam::PetscUtils::foamKSPMonitorRecordInit
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    void *cctx
)
{
    PetscFunctionBeginUser;
    if (it) PetscFunctionReturn(0);

    auto* ctx = static_cast<petscLinearSolverContext*>(cctx);
    solverPerformance& solverPerf = ctx->performance;

    solverPerf.initialResidual() = rnorm;
    PetscFunctionReturn(0);
}


PetscErrorCode Foam::PetscUtils::foamKSPConverge
(
    KSP ksp,
    PetscInt it,
    PetscReal rnorm,
    KSPConvergedReason* reason,
    void* cctx
)
{
    PetscFunctionBeginUser;
    //
    // Equivalent to the OpenFOAM checkConvergence function
    //
    auto* ctx = static_cast<petscLinearSolverContext*>(cctx);
    solverPerformance& solverPerf = ctx->performance;

    PetscReal rtol;
    PetscReal abstol;
    PetscReal divtol;
    PetscInt maxits;
    PetscCall(KSPGetTolerances(ksp, &rtol, &abstol, &divtol, &maxits));

    // compute L1 norm of residual (PETSc always uses L2)
    // assumes normFactor have been precomputed before solving the linear system
    // When using CG, this is actually a copy of the residual vector
    // stored inside the PETSc class (from PETSc version 3.14 on).
    // With GMRES instead, this call is more expensive
    // since we first need to generate the solution
    PetscCall(KSPBuildResidual(ksp, ctx->res_l1_w[0], ctx->res_l1_w[1], &ctx->res_l1_w[1]));
    PetscCall(VecNorm(ctx->res_l1_w[1], NORM_1, &rnorm));

    // rescale by the normFactor
    PetscReal residual = rnorm / ctx->normFactor;

    if (it == 0)
    {
        solverPerf.initialResidual() = residual;
    }

    if (residual < abstol)
    {
        *reason = KSP_CONVERGED_ATOL;
    }
    else if
    (
        rtol > Foam::SMALL  /** pTraits<Type>::one */
     && residual < rtol * solverPerf.initialResidual()
    )
    {
        *reason = KSP_CONVERGED_RTOL;
    }
    else if (it >= maxits)
    {
        *reason = KSP_DIVERGED_ITS;
    }
    else if (residual > divtol)
    {
        *reason = KSP_DIVERGED_DTOL;
    }
    solverPerf.finalResidual() = residual;

    PetscFunctionReturn(0);
}

void Foam::PetscUtils::setFlag
(
    const word& key,
    const word& val,
    const bool verbose
)
{
    if (verbose)
    {
        Info<< key << ' ' << val << nl;
    }

    AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), val.c_str()));
}


void Foam::PetscUtils::setFlags
(
    const word& prefix,
    const dictionary& dict,
    const bool verbose
)
{
    for (const entry& e : dict)
    {
        const word key = '-' + prefix + e.keyword();
        const word val = e.get<word>();

        if (verbose)
        {
            Info<< key << ' ' << val << nl;
        }

        AssertPETSc(PetscOptionsSetValue(NULL, key.c_str(), val.c_str()));
    }
}


// ************************************************************************* //
