/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "LduMatrix.H"
#include "diagTensorField.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
               *value;
        }
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solveSegregatedOrCoupled
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
            << "fvMatrix<Type>::solveSegregatedOrCoupled"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    // Do not solve if maxIter == 0
    if (solverControls.getOrDefault<label>("maxIter", -1) == 0)
    {
        return SolverPerformance<Type>();
    }

    word type(solverControls.getOrDefault<word>("type", "segregated"));

    if (type == "segregated")
    {
        return solveSegregated(solverControls);
    }
    else if (type == "coupled")
    {
        return solveCoupled(solverControls);
    }
    else
    {
        FatalIOErrorInFunction(solverControls)
            << "Unknown type " << type
            << "; currently supported solver types are segregated and coupled"
            << exit(FatalIOError);

        return SolverPerformance<Type>();
    }
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solveSegregated
(
    const dictionary& solverControls
)
{
    if (useImplicit_)
    {
        FatalErrorInFunction
            << "Implicit option is not allowed for type: " << Type::typeName
            << exit(FatalError);
    }

    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solveSegregated"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            SolverPerformance<Type>::debug
        );

    auto& psi =
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    SolverPerformance<Type> solverPerfVec
    (
        "fvMatrix<Type>::solveSegregated",
        psi.name()
    );

    scalarField saveDiag(diag());

    Field<Type> source(source_);

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implicit part by updateMatrixInterfaces within
    // the component loop.
    addBoundarySource(source);

    typename Type::labelType validComponents
    (
        psi.mesh().template validComponents<Type>()
    );

    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // copy field and source

        scalarField psiCmpt(psi.primitiveField().component(cmpt));
        addBoundaryDiag(diag(), cmpt);

        scalarField sourceCmpt(source.component(cmpt));

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        lduInterfaceFieldPtrsList interfaces =
            psi.boundaryField().scalarInterfaces();

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        {
            PrecisionAdaptor<solveScalar, scalar> sourceCmpt_ss(sourceCmpt);
            ConstPrecisionAdaptor<solveScalar, scalar> psiCmpt_ss(psiCmpt);

            const label startRequest = UPstream::nRequests();

            initMatrixInterfaces
            (
                true,
                bouCoeffsCmpt,
                interfaces,
                psiCmpt_ss(),
                sourceCmpt_ss.ref(),
                cmpt
            );

            updateMatrixInterfaces
            (
                true,
                bouCoeffsCmpt,
                interfaces,
                psiCmpt_ss(),
                sourceCmpt_ss.ref(),
                cmpt,
                startRequest
            );
        }

        solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        if (logLevel)
        {
            solverPerf.print(Info.masterStream(this->mesh().comm()));
        }

        solverPerfVec.replace(cmpt, solverPerf);
        solverPerfVec.solverName() = solverPerf.solverName();

        psi.primitiveFieldRef().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi.correctBoundaryConditions();

    psi.mesh().data().setSolverPerformance(psi.name(), solverPerfVec);

    return solverPerfVec;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solveCoupled
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info.masterStream(this->mesh().comm())
            << "fvMatrix<Type>::solveCoupled"
               "(const dictionary& solverControls) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            SolverPerformance<Type>::debug
        );

    auto& psi =
        const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

    LduMatrix<Type, scalar, scalar> coupledMatrix(psi.mesh());
    coupledMatrix.diag() = diag();
    coupledMatrix.upper() = upper();
    coupledMatrix.lower() = lower();
    coupledMatrix.source() = source();

    addBoundaryDiag(coupledMatrix.diag(), 0);
    addBoundarySource(coupledMatrix.source(), false);

    coupledMatrix.interfaces() = psi.boundaryFieldRef().interfaces();
    coupledMatrix.interfacesUpper() = boundaryCoeffs().component(0);
    coupledMatrix.interfacesLower() = internalCoeffs().component(0);

    autoPtr<typename LduMatrix<Type, scalar, scalar>::solver>
    coupledMatrixSolver
    (
        LduMatrix<Type, scalar, scalar>::solver::New
        (
            psi.name(),
            coupledMatrix,
            solverControls
        )
    );

    SolverPerformance<Type> solverPerf
    (
        coupledMatrixSolver->solve(psi)
    );

    if (logLevel)
    {
        solverPerf.print(Info.masterStream(this->mesh().comm()));
    }

    psi.correctBoundaryConditions();

    psi.mesh().data().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solveSegregatedOrCoupled()
{
    return this->solveSegregatedOrCoupled(solverDict());
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    return psi_.mesh().solve(*this, solverControls);
}


template<class Type>
Foam::autoPtr<typename Foam::fvMatrix<Type>::fvSolver>
Foam::fvMatrix<Type>::solver()
{
    return solver(solverDict());
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::fvSolver::solve()
{
    return solve(fvMat_.solverDict());
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solve(const word& name)
{
    return this->solve(solverDict(name));
}


template<class Type>
Foam::SolverPerformance<Type> Foam::fvMatrix<Type>::solve()
{
    return this->solve(solverDict());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::fvMatrix<Type>::residual() const
{
    auto tres = tmp<Field<Type>>::New(source_);
    auto& res = tres.ref();

    addBoundarySource(res);

    // Loop over field components
    for (direction cmpt=0; cmpt<Type::nComponents; cmpt++)
    {
        scalarField psiCmpt(psi_.primitiveField().component(cmpt));

        scalarField boundaryDiagCmpt(psi_.size(), Zero);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                psi_.boundaryField().scalarInterfaces(),
                cmpt
            )
        );
    }

    return tres;
}


// ************************************************************************* //
