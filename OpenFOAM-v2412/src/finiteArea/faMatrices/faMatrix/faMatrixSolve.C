/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

Description
    Finite-Area matrix basic solvers.

\*---------------------------------------------------------------------------*/

#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::faMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    internalCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

    boundaryCoeffs_[patchi][facei].component(cmpt) +=
        diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]*value;
}


template<class Type>
Foam::SolverPerformance<Type> Foam::faMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    DebugInFunction
        << "solving faMatrix<Type>"
        << endl;

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            SolverPerformance<Type>::debug
        );

    auto& psi =
        const_cast<GeometricField<Type, faPatchField, areaMesh>&>(psi_);

    SolverPerformance<Type> solverPerfVec
    (
        "faMatrix<Type>::solve",
        psi.name()
    );

    scalarField saveDiag(diag());

    Field<Type> source(source_);
    addBoundarySource(source);

    // Note: make a copy of interfaces: no longer a reference
    lduInterfaceFieldPtrsList interfaces =
        psi_.boundaryField().scalarInterfaces();

    for (direction cmpt = 0; cmpt < Type::nComponents; ++cmpt)
    {
        // copy field and source

        scalarField psiCmpt(psi_.primitiveField().component(cmpt));
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
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        if (logLevel)
        {
            solverPerf.print(Info);
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
Foam::SolverPerformance<Type> Foam::faMatrix<Type>::faSolver::solve()
{
    return solve(faMat_.solverDict());
}


template<class Type>
Foam::SolverPerformance<Type> Foam::faMatrix<Type>::solve(const word& name)
{
    return this->solve(solverDict(name));
}


template<class Type>
Foam::SolverPerformance<Type> Foam::faMatrix<Type>::solve()
{
    return this->solve(solverDict());
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::faMatrix<Type>::residual() const
{
    auto tres = tmp<Field<Type>>::New(source_);
    auto& res = tres().ref();

    addBoundarySource(res);

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; ++cmpt)
    {
        scalarField psiCmpt(psi_.internalField().component(cmpt));

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
