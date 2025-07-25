/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

Class
    Foam::petscLinearSolverContext

Description
    A struct containing a KSP and Mat for each equation to be solved

SourceFiles
    petscLinearSolverContext.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscLinearSolverContext_H
#define Foam_module_petscLinearSolverContext_H

#include "solverPerformance.H"
#include "petscCacheManager.H"
#include "petscErrorHandling.H"
#include "petscUtils.H"

// PETSc
#include "petscvec.h"
#include "petscmat.h"
#include "petscksp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class petscLinearSolverContext Declaration
\*---------------------------------------------------------------------------*/

class petscLinearSolverContext
{
    // Private Data

        bool init_;
        bool init_aux_vectors_;

public:

    // Public Members

        //- Using OpenFOAM norm for convergence testing
        bool useFoamTest;

        Mat Amat;
        KSP ksp;

        Vec sol, rhs;

        Vec ArowsSum;
        Vec *res_l1_w;

        solverPerformance performance;
        petscCacheManager caching;

        List<PetscInt> lowNonZero;
        PetscInt maxLowNonZeroPerRow = 0;

        PetscScalar normFactor;

        PetscLogStage matstage;
        PetscLogStage pcstage;
        PetscLogStage kspstage;


    // Constructors

        //- Default construct
        petscLinearSolverContext()
        :
            init_(false),
            init_aux_vectors_(false),
            useFoamTest(false),
            Amat(nullptr),
            ksp(nullptr),
            sol(nullptr),
            rhs(nullptr),
            ArowsSum(nullptr),
            res_l1_w(nullptr)
        {}


    //- Destructor
    virtual ~petscLinearSolverContext()
    {
        MatDestroy(&Amat);
        KSPDestroy(&ksp);
        VecDestroy(&sol);
        VecDestroy(&rhs);
        VecDestroy(&ArowsSum);
        if (res_l1_w) VecDestroyVecs(2,&res_l1_w);
    }


    // Member Functions

        //- Return value of initialized
        bool initialized() const noexcept
        {
            return init_;
        }

        //- Change state of initialized, return previous value
        bool initialized(const bool yes) noexcept
        {
            bool old(init_);
            init_ = yes;
            return old;
        }

        //- Create auxiliary rows for calculation purposes
        void createAuxVecs()
        {
            if (!Amat) return;
            if (!init_aux_vectors_)
            {
                AssertPETSc(MatCreateVecs(Amat, NULL, &ArowsSum));
                AssertPETSc(VecDuplicateVecs(ArowsSum, 2, &res_l1_w));
                init_aux_vectors_ = true;
            }
        }

        //- Sum the rows of A, placing the result in ArowsSum
        void initArowsSumVec()
        {
            if (init_aux_vectors_)
            {
                AssertPETSc(MatGetRowSum(Amat, ArowsSum));
            }
        }

        //- Compute the normFactor used in convergence testing,
        //- assumes ArowsSum has been already computed
        void computeNormFactor(Vec psi, Vec source)
        {
            if (!init_aux_vectors_) return;

            AssertPETSc(MatMult(Amat, psi, res_l1_w[1]));

            normFactor =
               Foam::PetscUtils::normFactor
               (
                   res_l1_w[1],
                   psi,
                   source,
                   ArowsSum
               );
        }

        //- Compute residual (L1) norm
        //- assumes normFactor has been already computed
        PetscReal getResidualNormL1(Vec psi, Vec source)
        {
            if (!init_aux_vectors_) return static_cast<PetscReal>(-1);

            PetscReal rnorm;
            AssertPETSc(MatMult(Amat, psi, res_l1_w[1]));
            AssertPETSc(VecAXPY(res_l1_w[1], -1., source));
            AssertPETSc(VecNorm(res_l1_w[1], NORM_1, &rnorm));
            return rnorm/normFactor;
        }

        //- Compute residual (L1) norm
        //- assumes normFactor has been already computed
        //- as well as Adotpsi is stored in res_l1_w[1]
        PetscReal getResidualNormL1(Vec source)
        {
            if (!init_aux_vectors_) return static_cast<PetscReal>(-1);

            PetscReal rnorm;
            AssertPETSc(VecAXPY(res_l1_w[1], -1., source));
            AssertPETSc(VecNorm(res_l1_w[1], NORM_1, &rnorm));
            return rnorm/normFactor;
        }


    // Housekeeping

        //- Non-const access to state of initialized
        bool& initialized() noexcept
        {
            return init_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
