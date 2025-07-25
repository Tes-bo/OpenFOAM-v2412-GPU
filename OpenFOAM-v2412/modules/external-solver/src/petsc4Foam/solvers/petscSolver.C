/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
    Copyright (C) 2019-2020 Simone Bna
    Copyright (C) 2021 Stefano Zampini
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

#include "fvMesh.H"
#include "fvMatrices.H"
#include "globalIndex.H"
#include "PrecisionAdaptor.H"
#include "cyclicLduInterface.H"
#include "cyclicAMILduInterface.H"
#include "PstreamGlobals.H"
#include "addToRunTimeSelectionTable.H"

#include "petscSolver.H"
#include "petscControls.H"
#include "petscLinearSolverContexts.H"
#include "petscUtils.H"
#include "petscErrorHandling.H"
#include "petscWrappedVector.H"

#include <algorithm>  // For std::max_element
#include <cstring>    // For NULL

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(petscSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<petscSolver>
        addpetscSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<petscSolver>
        addpetscSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

// Test if given matrix appears to be an "sbaij" type
static bool is_sbaij(Mat& Amat)
{
    MatType mtype;
    AssertPETSc(MatGetType(Amat, &mtype));

    return
    (
        !strcmp(mtype, "sbaij")
     || !strcmp(mtype, "seqsbaij")
     || !strcmp(mtype, "mpisbaij")
    );
}

} // End anonymous namespace

#if PETSC_VERSION_LT(3,17,0)
typedef PetscInt64 PetscCount;
#endif

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::petscSolver::petscSolver
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
    ),
    petscDict_(solverControls.subDict("petsc")),
    eqName_(fieldName),
    prefix_("eqn_" + eqName_ + "_")
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::solverPerformance Foam::petscSolver::scalarSolve
(
    solveScalarField& psi,
    const solveScalarField& source,
    const direction cmpt
) const
{
    // Ensure PETSc is initialized
    const petscControls& controls = petscControls::New
    (
        matrix_.mesh().thisDb().time()
    );

    if (!controls.valid())
    {
        FatalErrorInFunction
            << "PETSc not initialized" << nl << endl;
    }

    const bool verbose
    (
        petscDict_.getOrDefault("verbose", false)
    );

    // set all petsc flags in the prefix db
    dictionary petscDictOptions = petscDict_.subOrEmptyDict("options");
    PetscUtils::setFlags
    (
        prefix_,
        petscDictOptions,
        verbose
    );

    const fvMesh& fvm = dynamicCast<const fvMesh>(matrix_.mesh().thisDb());

    const petscLinearSolverContexts& contexts =
        petscLinearSolverContexts::New(fvm);

    petscLinearSolverContext& ctx = contexts.getContext(eqName_);
    const bool firsttimein = !ctx.initialized();

    if (firsttimein)
    {
        dictionary petscDictCaching = petscDict_.subOrEmptyDict("caching");
        ctx.caching.init(petscDictCaching);
    }
    ctx.caching.eventBegin();

    Mat& Amat = ctx.Amat;
    KSP& ksp  = ctx.ksp;
    Vec& rhs  = ctx.rhs;
    Vec& sol  = ctx.sol;

    List<PetscInt>& lowNonZero = ctx.lowNonZero;
    PetscInt& maxLowNonZeroPerRow = ctx.maxLowNonZeroPerRow;

    if (firsttimein)
    {
        DebugInfo<< "Initializing PETSc Linear Solver " << eqName_ << nl;

        ctx.initialized(true);

        AssertPETSc(PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_mat").c_str(),
            &ctx.matstage
        ));
        AssertPETSc(PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_pc").c_str(),
            &ctx.pcstage
        ));
        AssertPETSc(PetscLogStageRegister
        (
            ("foam_" + eqName_ + "_ksp").c_str(),
            &ctx.kspstage
        ));
        AssertPETSc(PetscLogStagePush(ctx.matstage));
        buildMat(Amat, lowNonZero, maxLowNonZeroPerRow);
        buildKsp(Amat, ksp);
        // Create solution and rhs vectors for PETSc only once
        AssertPETSc(MatCreateVecs(Amat, &sol, &rhs));
        AssertPETSc(PetscLogStagePop());
    }

    const bool matup = ctx.caching.needsMatrixUpdate();
    const bool pcup = ctx.caching.needsPrecondUpdate();
    if (matup || firsttimein)
    {
        DebugInfo<< "Update PETSc Linear Solver Matrix for: " << eqName_ << nl;
        AssertPETSc(PetscLogStagePush(ctx.matstage));
        updateMat(Amat, lowNonZero, maxLowNonZeroPerRow);
        AssertPETSc(PetscLogStagePop());
    }
    updateKsp(ksp, Amat, !pcup);

    // Use built-in residual norm computation
    const bool usePetscResidualNorm
    (
        petscDict_.getOrDefault("use_petsc_residual_norm", false)
    );

    // Monitor built-in residual
    const bool monitorFoamResidualNorm
    (
        petscDict_.getOrDefault("monitor_foam_residual_norm", false)
    );

    // This optimization is disabled since some KSP implementations
    // are buggy in PETSc wrt KSP_NORM_NONE
    // users can still provide -ksp_norm_type none at command line
#if 0
    // Disable KSP default computation of residual norm if we are
    // monitoring convergence a-la OpenFOAM
    if (!usePetscResidualNorm)
    {
        KSPSetNormType(ksp, KSP_NORM_NONE);
    }
    else
    {
        KSPSetNormType(ksp, KSP_NORM_DEFAULT);
    }
#endif

    // ksp set options from db (may change norm type here if needed)
    if (firsttimein) AssertPETSc(KSPSetFromOptions(ksp));

    // Solver name from petsc (may have been changed from option database)
    KSPType ksptype;
    AssertPETSc(KSPGetType(ksp, &ksptype));
    const word solverName(ksptype);

    // Setup class containing solver performance data
    solverPerformance solverPerf
    (
        "PETSc-" + solverName,
        fieldName_
    );

    // Retain copy of solverPerformance
    ctx.performance = solverPerf;

    // Place openfoam raw memory in PETSc vectors
    PetscWrappedVector petsc_psi(psi, sol);
    PetscWrappedVector petsc_source(source, rhs);

    // Initialize data to compute L1 residual or to monitor it while solving
    ctx.performance.initialResidual() = 0;
    KSPNormType normType;
    AssertPETSc(KSPGetNormType(ksp, &normType));
    if (!usePetscResidualNorm || monitorFoamResidualNorm || normType == KSP_NORM_NONE)
    {
        ctx.createAuxVecs();
        ctx.initArowsSumVec();
        ctx.computeNormFactor(petsc_psi, petsc_source);
        if (normType == KSP_NORM_NONE && !usePetscResidualNorm)
        {
            ctx.performance.initialResidual() =
                ctx.getResidualNormL1(petsc_source);
        }
    }

    // Convergence testing
    if (firsttimein)
    {
        if (usePetscResidualNorm)
        {
            // Add monitor to record initial residual computed by PETSc
            AssertPETSc(KSPMonitorSet
            (
                ksp,
                PetscUtils::foamKSPMonitorRecordInit,
                &ctx,
                NULL
            ));
        }
        else
        {
            ctx.useFoamTest = true;
            AssertPETSc(KSPSetConvergenceTest
            (
                ksp,
                PetscUtils::foamKSPConverge,
                &ctx,
                NULL
            ));
        }

        // Monitoring (print to stdout) residual reduction in OpenFOAM norm
        if (monitorFoamResidualNorm)
        {
            AssertPETSc(KSPMonitorSet
            (
                ksp,
                PetscUtils::foamKSPMonitorFoam,
                &ctx,
                NULL
            ));
        }
    }

    // Geometric multigrid support
    // (pc_type ml or pc_type gamg pc_gamg_type geo for 2d only as for 3.14)
    if (firsttimein)
    {
        PC pc;
        void (*f)(void) = NULL;
        PetscBool hypre_device = PETSC_FALSE;

        #ifdef PETSC_HAVE_HYPRE_DEVICE
        AssertPETSc(PetscObjectTypeCompare((PetscObject)pc, PCHYPRE, &hypre_device));
        #endif
        AssertPETSc(KSPGetPC(ksp, &pc));
        PetscObjectQueryFunction((PetscObject)pc, "PCSetCoordinates_C", &f);

        // Cell-centres are contiguous in memory (array of structures)
        // Not sure about const-ness here
        // Could also use the PrecisionAdaptor
        if (f && !hypre_device)
        {
            const vectorField& cc = fvm.C().primitiveField();
            const PetscInt n = cc.size();
            const PetscInt sdim = vector::nComponents;

            List<PetscReal> ccPoints(cc.size()*vector::nComponents);

            auto iter = ccPoints.data();
            for (const vector& v : cc)
            {
                *(iter++) = v.x();
                *(iter++) = v.y();
                *(iter++) = v.z();
            }

            AssertPETSc(PCSetCoordinates(pc, sdim, n, ccPoints.data()));
        }
    }

    // Setup KSP (this is not explicitly needed, but we call it to separate
    // PCSetUp from KSPSolve timings when requesting -log_view from PETSc)
    AssertPETSc(PetscLogStagePush(ctx.pcstage));
    AssertPETSc(KSPSetUp(ksp));
    AssertPETSc(KSPSetUpOnBlocks(ksp));
    AssertPETSc(PetscLogStagePop());

    // Fail if setup failed
    KSPConvergedReason reason;
    AssertPETSc(KSPGetConvergedReason(ksp, &reason));
    if (reason)
    {
       FatalErrorInFunction
           << "Linear solver setup failed: reason "
           << KSPConvergedReasons[reason]
           << nl
           << exit(FatalError);
    }

    // Solve A x = b
    AssertPETSc(PetscLogStagePush(ctx.kspstage));
    AssertPETSc(KSPSolve(ksp, petsc_source, petsc_psi));
    AssertPETSc(PetscLogStagePop());

    ctx.caching.eventEnd();

    // Set nIterations and final residual
    PetscInt nIters;
    AssertPETSc(KSPGetIterationNumber(ksp, &nIters));
    ctx.performance.nIterations() = nIters;

    // If using the PETSc infrastructure to monitor convergence,
    // we need to report back the final residual
    if (usePetscResidualNorm)
    {
        PetscReal rnorm;
        KSPNormType normType;
        AssertPETSc(KSPGetNormType(ksp, &normType));
        if (normType == KSP_NORM_NONE)
        {
            // We run with PETSc norm KSP_NORM_NONE -> report L1 norm
            rnorm = ctx.getResidualNormL1(petsc_psi, petsc_source);
        }
        else
        {
            // Report final PETSc norm since users explicitly ask to run
            // with PETSc norm
            AssertPETSc(KSPGetResidualNorm(ksp, &rnorm));
        }
        ctx.performance.finalResidual() = rnorm;
    }

    // Return solver performance to OpenFOAM
    return ctx.performance;
}


Foam::solverPerformance Foam::petscSolver::solve
(
    scalarField& psi_s,
    const scalarField& source,
    const direction cmpt
) const
{
    PrecisionAdaptor<solveScalar, scalar> tpsi(psi_s);
    return scalarSolve
    (
        tpsi.ref(),
        ConstPrecisionAdaptor<solveScalar, scalar>(source)(),
        cmpt
    );
}


void Foam::petscSolver::updateKsp
(
    KSP& ksp,
    Mat& Amat,
    const bool reuse
) const
{
    if (reuse)
    {
        DebugInfo<< "Cache-Hit: reuse preconditioner " << eqName_ << nl;

        AssertPETSc(KSPSetReusePreconditioner(ksp, PETSC_TRUE));
    }
    else
    {
        DebugInfo<< "Cache-Hit: rebuild preconditioner " << eqName_ << nl;

        AssertPETSc(KSPSetReusePreconditioner(ksp, PETSC_FALSE));
    }

    AssertPETSc(KSPSetOperators(ksp, Amat, Amat));

    // update tolerance and relTol for the (*)Final solver
    AssertPETSc(KSPSetTolerances
    (
        ksp,
        relTol_,
        tolerance_,
        PETSC_DEFAULT,
        maxIter_
    ));
}


void Foam::petscSolver::buildKsp
(
    Mat& Amat,
    KSP& ksp
) const
{
    // Create parallel solver context
    AssertPETSc(KSPCreate(PetscObjectComm((PetscObject)Amat), &ksp));

    // Set the prefix for the options db (e.g. -eqn_p_)
    AssertPETSc(KSPSetOptionsPrefix(ksp, prefix_.c_str()));

    // ksp set operator and preconditioner
    AssertPETSc(KSPSetOperators(ksp, Amat, Amat));

    // OpenFOAM relative tolerance -> tolerance_
    AssertPETSc(KSPSetTolerances
    (
        ksp,
        relTol_,
        tolerance_,
        PETSC_DEFAULT,
        maxIter_
    ));

    // Use solution from the previous timestep as initial guess
    AssertPETSc(KSPSetInitialGuessNonzero(ksp, PETSC_TRUE));
}


void Foam::petscSolver::buildMat
(
    Mat& Amat,
    List<PetscInt>& lowNonZero,
    PetscInt& maxLowNonZeroPerRow
) const
{
    // communicator for the matrix
    // PETSc internally duplicates the communicator to
    // avoid tag/attributes clashing

    const auto comm =
    (
        PstreamGlobals::MPICommunicators_.empty()
      ? PETSC_COMM_SELF
      : PstreamGlobals::MPICommunicators_[matrix_.mesh().comm()]
    );

    // traversal-based algorithm
    typedef List<PetscInt> integerList;

    const lduAddressing& lduAddr = matrix_.mesh().lduAddr();
    const lduInterfacePtrsList interfaces(matrix_.mesh().interfaces());

    const labelUList& upp = lduAddr.upperAddr();
    const labelUList& low = lduAddr.lowerAddr();
    // const labelUList& losort = lduAddr.losortAddr();

    // Local degrees-of-freedom i.e. number of local rows
    const label nrows_ = lduAddr.size();

    // Create a petsc matrix
    AssertPETSc(MatCreate(comm, &Amat));
    AssertPETSc(MatSetSizes(Amat, nrows_, nrows_, PETSC_DECIDE, PETSC_DECIDE));

    // Set the prefix for the options db (e.g. -eqn_p_)
    AssertPETSc(MatSetOptionsPrefix(Amat, prefix_.c_str()));
    AssertPETSc(MatSetFromOptions(Amat));

    // Query for fast COO support
    #if PETSC_VERSION_LT(3,17,0)
    PetscErrorCode (*f)(Mat,PetscInt,const PetscInt[],const PetscInt[]) = NULL;
    #else
    PetscErrorCode (*f)(Mat,PetscCount,const PetscInt[],const PetscInt[]) = NULL;
    #endif
    AssertPETSc(PetscObjectQueryFunction((PetscObject)Amat,"MatSetPreallocationCOO_C",&f));
    if (f)
    {
        const globalIndex globalNumbering_(nrows_);

        // Number of non-zeros from interfaces (on-processor or off-processor)
        label nProcValues = 0;
        labelList globalCells
        (
            identity(globalNumbering_.range())
        );

        const label startOfRequests = UPstream::nRequests();
        forAll(interfaces, patchi)
        {
            const lduInterface* intf = interfaces.set(patchi);

            if (intf)
            {
                if (isA<cyclicAMILduInterface>(*intf))
                {
                    FatalErrorInFunction
                        << "cyclicAMI is not supported" << nl
                        << exit(FatalError);
                }
                nProcValues += lduAddr.patchAddr(patchi).size();
                interfaces[patchi].initInternalFieldTransfer
                (
                    Pstream::commsTypes::nonBlocking,
                    globalCells
                );
            }
        }
        PetscCount toti = nrows_+low.size()+upp.size()+nProcValues;

        PetscInt *coo_i, *coo_j;
        AssertPETSc(PetscMalloc1(toti, &coo_i));
        AssertPETSc(PetscMalloc1(toti, &coo_j));

        PetscInt *ptr_i = coo_i;
        PetscInt *ptr_j = coo_j;
        label off = globalNumbering_.localStart();
        for (label celli = off; celli < off + nrows_; ++celli)
        {
            const PetscInt globRow = celli;
            *ptr_i++ = globRow;
            *ptr_j++ = globRow;
        }
        for (label idx=0; idx < upp.size(); ++idx)
        {
            PetscInt globRow = low[idx] + off;
            PetscInt globCol = upp[idx] + off;
            *ptr_i++ = globRow;
            *ptr_j++ = globCol;
        }
        for (label idx=0; idx < low.size(); ++idx)
        {
            PetscInt globRow = upp[idx] + off;
            PetscInt globCol = low[idx] + off;
            *ptr_i++ = globRow;
            *ptr_j++ = globCol;
        }
        if (UPstream::parRun())
        {
            UPstream::waitRequests(startOfRequests);
        }
        forAll(interfaces, patchi)
        {
            if (interfaces.set(patchi))
            {
                // Processor-local values
                const labelUList& faceCells = lduAddr.patchAddr(patchi);
                labelField nbrCells
                (
                    interfaces[patchi].internalFieldTransfer
                    (
                        Pstream::commsTypes::nonBlocking,
                        globalCells
                    )
                );

                if (faceCells.size() != nbrCells.size())
                {
                    FatalErrorInFunction
                        << "Mismatch in interface sizes (AMI?)" << nl
                        << "Have " << faceCells.size() << " != "
                        << nbrCells.size() << nl
                        << exit(FatalError);
                }

                forAll(faceCells, i)
                {
                    PetscInt globRow = faceCells[i] + off;
                    PetscInt globCol = nbrCells[i];
                    *ptr_i++ = globRow;
                    *ptr_j++ = globCol;
                }
            }
        }

        // Preallocate with COO information
        // Use interface call to get the timing logged
#if PETSC_VERSION_GE(3,15,0)
        AssertPETSc(MatSetPreallocationCOO(Amat,toti,coo_i,coo_j));
#else
        (*f)(Amat,toti,coo_i,coo_j);
#endif
        AssertPETSc(PetscFree(coo_i));
        AssertPETSc(PetscFree(coo_j));
    }
    else
    {
        // On-processor non-zero entries, per row.
        // Initialise with 1 (includes diagonal).
        integerList nonZero_(nrows_, 1);
        lowNonZero.resize(nrows_, 0);

        // Number of on-processor non-zeros
        // set the matrix options

        if (!is_sbaij(Amat))
        {
            for (const label celli : upp)
            {
                ++nonZero_[celli];
            }
        }

        for (const label celli : low)
        {
            ++nonZero_[celli];
            ++lowNonZero[celli];
        }


        // PetscInt may be an unusually sized type.
        // Use max_element instead of Foam::max(List..)

        maxLowNonZeroPerRow = 0;
        {
            auto maxElem =
                std::max_element(lowNonZero.cbegin(), lowNonZero.cend());

            if (maxElem != lowNonZero.cend())  // Extra safety
            {
                maxLowNonZeroPerRow = *maxElem;
            }
        }

        // Off-processor non-zero entries, per row.

        integerList nonZeroNeigh_(nrows_, 0);

        // Number of non-zeros from interfaces (on-processor or off-processor)
        forAll(interfaces, patchi)
        {
            const lduInterface* intf = interfaces.set(patchi);

            if (intf)
            {
                const labelUList& faceCells = lduAddr.patchAddr(patchi);

                if (isA<cyclicAMILduInterface>(*intf))
                {
                    FatalErrorInFunction
                        << "cyclicAMI is not supported" << nl
                        << exit(FatalError);
                }
                else if (isA<cyclicLduInterface>(*intf))
                {
                    // Cyclic is on-processor
                    for (const label celli : faceCells)
                    {
                        ++nonZero_[celli];
                    }
                }
                else
                {
                    // Off-processor
                    for (const label celli : faceCells)
                    {
                        ++nonZeroNeigh_[celli];
                    }
                }
            }
        }

        // preallocate the matrix
        AssertPETSc(MatXAIJSetPreallocation
        (
            Amat,
            1,
            nonZero_.data(),        // on-processor
            nonZeroNeigh_.data(),   // off-processor
            nonZero_.data(),
            nonZeroNeigh_.data()
        ));

        // The general AIJ preallocator does not handle SELL
        AssertPETSc(MatSeqSELLSetPreallocation
        (
            Amat,
            0, nonZero_.data()        // on-processor
        ));

        AssertPETSc(MatMPISELLSetPreallocation
        (
            Amat,
            0, nonZero_.data(),       // on-processor
            0, nonZeroNeigh_.data()   // off-processor
        ));
    }

    // set the matrix options
    if (matrix_.symmetric())
    {
        AssertPETSc(MatSetOption(Amat, MAT_SYMMETRIC, PETSC_TRUE));
        AssertPETSc(MatSetOption(Amat, MAT_SYMMETRY_ETERNAL, PETSC_TRUE));
    }
    else
    {
        AssertPETSc(MatSetOption(Amat, MAT_STRUCTURALLY_SYMMETRIC, PETSC_TRUE));
    }

    AssertPETSc(MatSetOption(Amat, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE));

    // Skip MPI_Allreduce calls when finalizing assembly
    AssertPETSc(MatSetOption(Amat, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE));

    // Allow general setup of other matrices instead of erroring
    AssertPETSc(MatSetUp(Amat));
}


void Foam::petscSolver::updateMat
(
    Mat& Amat,
    const List<PetscInt>& lowNonZero,
    const PetscInt maxLowNonZeroPerRow
) const
{
    const lduAddressing& lduAddr = matrix_.mesh().lduAddr();
    const lduInterfacePtrsList interfaces(matrix_.mesh().interfaces());

    const labelUList& upp = lduAddr.upperAddr();
    const labelUList& low = lduAddr.lowerAddr();
    // const labelUList& losort = lduAddr.losortAddr();

    const scalarField& diagVal = matrix_.diag();
    const scalarField& uppVal = matrix_.upper();
    const scalarField& lowVal =
    (
        matrix_.hasLower()
      ? matrix_.lower()
      : matrix_.upper()
    );

    // Local degrees-of-freedom i.e. number of local rows
    const label nrows_ = lduAddr.size();

    // Query for fast COO support
    PetscErrorCode (*f)(Mat,const PetscScalar[],InsertMode) = NULL;
    AssertPETSc(PetscObjectQueryFunction((PetscObject)Amat,"MatSetValuesCOO_C",&f));
    if (f)
    {
       label nProcValues = 0;
       forAll(interfaces, patchi)
       {
          const lduInterface* intf = interfaces.set(patchi);

          if (intf)
          {
             nProcValues += lduAddr.patchAddr(patchi).size();
          }
       }
       PetscCount toti = nrows_+low.size()+upp.size()+nProcValues;

       PetscScalar *coo_v;
       PetscMalloc1(toti,&coo_v);
       PetscArraycpy(coo_v,diagVal.cdata(),nrows_);
       PetscArraycpy(coo_v+nrows_,uppVal.cdata(),upp.size());
       PetscArraycpy(coo_v+nrows_+upp.size(),lowVal.cdata(),low.size());
       if (nProcValues)
       {
          PetscScalar *ptr_v = coo_v+nrows_+low.size()+upp.size();
          forAll(interfaces, patchi)
          {
             const lduInterface* intf = interfaces.set(patchi);

             if (intf)
             {
                const labelUList& faceCells = lduAddr.patchAddr(patchi);
                const scalarField& bCoeffs = interfaceBouCoeffs_[patchi];
                // NB:opposite sign since we using this sides'
                //    coefficients (from discretising our face; not the
                //    neighbouring, reversed face)
                forAll(faceCells, i)
                {
                   *ptr_v++ = -bCoeffs[i];
                }
             }
          }

       }
#if PETSC_VERSION_GE(3,15,0)
       AssertPETSc(MatSetValuesCOO(Amat,coo_v,INSERT_VALUES));
#else
       (*f)(Amat,coo_v,INSERT_VALUES);
#endif
       AssertPETSc(PetscFree(coo_v));
    }
    else
    {
       // Number of internal faces (connectivity)
       const label nIntFaces_ = upp.size();

       const globalIndex globalNumbering_(nrows_);

       // The diagonal
       {
           PetscInt globRow = globalNumbering_.toGlobal(0);

           for (label celli = 0; celli < nrows_; ++celli)
           {
               PetscScalar val = diagVal[celli];

               AssertPETSc(MatSetValues
               (
                   Amat,
                   1, &globRow, // row
                   1, &globRow, // col
                   &val,
                   INSERT_VALUES
               ));
               ++globRow;
           }
       }

       // Connections to higher numbered cells (on-processor) (optimized form)
       PetscInt globRow;
       PetscInt globCols[maxLowNonZeroPerRow];
       PetscScalar globVals[maxLowNonZeroPerRow];
       label kidx = 0;

       for (label idx=0; idx < nrows_; ++idx)
       {
           if (lowNonZero[idx] > 0)
           {
               globRow = globalNumbering_.toGlobal(low[kidx]);

               for (label jidx=0; jidx < lowNonZero[idx]; ++jidx)
               {
                   // The row is sorted, col is unsorted
                   globCols[jidx] = globalNumbering_.toGlobal(upp[kidx]);
                   globVals[jidx] = uppVal[kidx];
                   ++kidx;
               }

               AssertPETSc(MatSetValues
               (
                   Amat,
                   1, &globRow,
                   lowNonZero[idx], globCols,
                   globVals,
                   INSERT_VALUES
               ));
           }
       }

       if (!is_sbaij(Amat))
       {
           // Connections to lower numbered cells (on-processor)
           for (label idx=0; idx < nIntFaces_; ++idx)
           {
               // The row is unsorted, col is sorted
               PetscInt globRow = globalNumbering_.toGlobal(upp[idx]);
               PetscInt globCol = globalNumbering_.toGlobal(low[idx]);
               PetscScalar val = lowVal[idx];

               AssertPETSc(MatSetValues
               (
                   Amat,
                   1, &globRow,
                   1, &globCol,
                   &val,
                   INSERT_VALUES
               ));
           }
       }

       labelList globalCells
       (
           identity
           (
               globalNumbering_.localSize(),
               globalNumbering_.localStart()
           )
       );

       // Connections to neighbouring processors
       {
           const label startOfRequests = UPstream::nRequests();

           // Initialise transfer of global cells
           forAll(interfaces, patchi)
           {
               if (interfaces.set(patchi))
               {
                   interfaces[patchi].initInternalFieldTransfer
                   (
                       Pstream::commsTypes::nonBlocking,
                       globalCells
                   );
               }
           }

           if (UPstream::parRun())
           {
               UPstream::waitRequests(startOfRequests);
           }

           forAll(interfaces, patchi)
           {
               if (interfaces.set(patchi))
               {
                   // Processor-local values
                   const labelUList& faceCells = lduAddr.patchAddr(patchi);

                   const scalarField& bCoeffs = interfaceBouCoeffs_[patchi];

                   labelField nbrCells
                   (
                       interfaces[patchi].internalFieldTransfer
                       (
                           Pstream::commsTypes::nonBlocking,
                           globalCells
                       )
                   );


                   if (faceCells.size() != nbrCells.size())
                   {
                       FatalErrorInFunction
                           << "Mismatch in interface sizes (AMI?)" << nl
                           << "Have " << faceCells.size() << " != "
                           << nbrCells.size() << nl
                           << exit(FatalError);
                   }

                   const label off = globalNumbering_.localStart();

                   // NB:opposite sign since we using this sides'
                   //    coefficients (from discretising our face; not the
                   //    neighbouring, reversed face)

                   forAll(faceCells, i)
                   {
                       PetscInt globRow = faceCells[i] + off; // row is sorted
                       PetscInt globCol = nbrCells[i];  // <- already global
                       PetscScalar val = -bCoeffs[i];

                       AssertPETSc(MatSetValues
                       (
                           Amat,
                           1, &globRow,
                           1, &globCol,
                           &(val),
                           INSERT_VALUES
                       ));
                   }
               }
           }
       }
       AssertPETSc(MatAssemblyBegin(Amat, MAT_FINAL_ASSEMBLY));
       AssertPETSc(MatAssemblyEnd(Amat, MAT_FINAL_ASSEMBLY));
    }
}


// ************************************************************************* //
