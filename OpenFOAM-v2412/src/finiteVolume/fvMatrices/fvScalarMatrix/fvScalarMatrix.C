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

#include "fvScalarMatrix.H"
#include "extrapolatedCalculatedFvPatchFields.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"
#include "jumpCyclicFvPatchField.H"
#include "cyclicPolyPatch.H"
#include "cyclicAMIPolyPatch.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvMatrix<Foam::scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
               *value;
        }
    }
}


template<>
Foam::autoPtr<Foam::fvMatrix<Foam::scalar>::fvSolver>
Foam::fvMatrix<Foam::scalar>::solver
(
    const dictionary& solverControls
)
{
    word regionName;
    if (psi_.mesh().name() != polyMesh::defaultRegion)
    {
        regionName = psi_.mesh().name() + "::";
    }
    addProfiling(solve, "fvMatrix::solve.", regionName, psi_.name());

    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solver(const dictionary& solverControls) : "
               "solver for fvMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    lduInterfaceFieldPtrsList interfaces =
        psi_.boundaryField().scalarInterfaces();


    autoPtr<fvMatrix<scalar>::fvSolver> solverPtr
    (
        new fvMatrix<scalar>::fvSolver
        (
            *this,
            lduMatrix::solver::New
            (
                psi_.name(),
                *this,
                boundaryCoeffs_,
                internalCoeffs_,
                interfaces,
                solverControls
            )
        )
    );

    diag() = saveDiag;

    return solverPtr;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::fvSolver::solve
(
    const dictionary& solverControls
)
{
    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            solverPerformance::debug
        );

    auto& psi =
        const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
        (
            fvMat_.psi()
        );

    scalarField saveDiag(fvMat_.diag());
    fvMat_.addBoundaryDiag(fvMat_.diag(), 0);

    scalarField totalSource(fvMat_.source());
    fvMat_.addBoundarySource(totalSource, false);

    // Assign new solver controls
    solver_->read(solverControls);

    solverPerformance solverPerf = solver_->solve
    (
        psi.primitiveFieldRef(),
        totalSource
    );

    if (logLevel)
    {
        solverPerf.print(Info.masterStream(fvMat_.mesh().comm()));
    }

    fvMat_.diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().data().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::solverPerformance Foam::fvMatrix<Foam::scalar>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<scalar>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            solverPerformance::debug
        );

    scalarField saveLower;
    scalarField saveUpper;

    if (useImplicit_)
    {
        createOrUpdateLduPrimitiveAssembly();

        if (psi_.mesh().fluxRequired(psi_.name()))
        {
            // Save lower/upper for flux calculation
            if (asymmetric())
            {
                saveLower = lower();
            }
            saveUpper = upper();
        }

        setLduMesh(*lduMeshPtr());
        transferFvMatrixCoeffs();
        setBounAndInterCoeffs();
        direction cmpt = 0;
        manipulateMatrix(cmpt);
    }

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    scalarField totalSource(source_);
    addBoundarySource(totalSource, false);

    lduInterfaceFieldPtrsList interfaces;
    PtrDynList<lduInterfaceField> newInterfaces;
    if (!useImplicit_)
    {
        interfaces = this->psi(0).boundaryField().scalarInterfaces();
    }
    else
    {
        setInterfaces(interfaces, newInterfaces);
    }

    tmp<scalarField> tpsi;
    if (!useImplicit_)
    {
        tpsi.ref
        (
            const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
            (
                psi_
            ).primitiveFieldRef()
        );
    }
    else
    {
        tpsi = tmp<scalarField>::New(lduAddr().size(), Zero);
        scalarField& psi = tpsi.ref();

        for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
        {
            const label cellOffset = lduMeshPtr()->cellOffsets()[fieldi];
            const auto& psiInternal = this->psi(fieldi).primitiveField();

            forAll(psiInternal, localCellI)
            {
                psi[cellOffset + localCellI] = psiInternal[localCellI];
            }
        }
    }
    scalarField& psi = tpsi.ref();

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        this->psi(0).name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces,
        solverControls
    )->solve(psi, totalSource);

    if (useImplicit_)
    {
        for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
        {
            auto& psiInternal =
                const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
                (
                    this->psi(fieldi)
                ).primitiveFieldRef();

            const label cellOffset = lduMeshPtr()->cellOffsets()[fieldi];

            forAll(psiInternal, localCellI)
            {
                psiInternal[localCellI] = psi[localCellI + cellOffset];
            }
        }
    }

    if (logLevel)
    {
        solverPerf.print(Info.masterStream(mesh().comm()));
    }

    diag() = saveDiag;

    if (useImplicit_)
    {
        if (psi_.mesh().fluxRequired(psi_.name()))
        {
            // Restore lower/upper
            if (asymmetric())
            {
                lower().setSize(saveLower.size());
                lower() = saveLower;
            }

            upper().setSize(saveUpper.size());
            upper() = saveUpper;
        }
        // Set the original lduMesh
        setLduMesh(psi_.mesh());
    }

    for (label fieldi = 0; fieldi < nMatrices(); fieldi++)
    {
        auto& localPsi =
            const_cast<GeometricField<scalar, fvPatchField, volMesh>&>
            (
                this->psi(fieldi)
            );

        localPsi.correctBoundaryConditions();
        localPsi.mesh().data().setSolverPerformance
        (
            localPsi.name(),
            solverPerf
        );
    }

    return solverPerf;
}


template<>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Foam::scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), Zero);
    addBoundaryDiag(boundaryDiag, 0);

    const scalarField& psif = psi_.primitiveField();
    ConstPrecisionAdaptor<solveScalar, scalar> tpsi(psif);
    const solveScalarField& psi = tpsi();

    tmp<solveScalarField> tres
    (
        lduMatrix::residual
        (
            psi,
            source_ - boundaryDiag*psif,
            boundaryCoeffs_,
            psi_.boundaryField().scalarInterfaces(),
            0
        )
    );

    ConstPrecisionAdaptor<scalar, solveScalar> tres_s(tres);
    addBoundarySource(tres_s.ref());

    return tres_s;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H() const
{
    auto tHphi = volScalarField::New
    (
        "H(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/dimVol,
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Hphi = tHphi.ref();

    Hphi.primitiveFieldRef() = (lduMatrix::H(psi_.primitiveField()) + source_);
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.primitiveFieldRef() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H1() const
{
    auto tH1 = volScalarField::New
    (
        "H(1)",
        psi_.mesh(),
        dimensions_/(dimVol*psi_.dimensions()),
        fvPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& H1_ = tH1.ref();

    H1_.primitiveFieldRef() = lduMatrix::H1();
    //addBoundarySource(Hphi.primitiveField());

    H1_.primitiveFieldRef() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


// ************************************************************************* //
