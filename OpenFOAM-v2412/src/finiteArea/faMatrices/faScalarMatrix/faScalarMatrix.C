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

\*---------------------------------------------------------------------------*/

#include "faScalarMatrix.H"
#include "extrapolatedCalculatedFaPatchFields.H"
#include "profiling.H"
#include "PrecisionAdaptor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::faMatrix<Foam::scalar>::setComponentReference
(
    const label patchI,
    const label edgeI,
    const direction,
    const scalar value
)
{
    const labelUList& faceLabels = psi_.mesh().boundary()[patchI].edgeFaces();

    internalCoeffs_[patchI][edgeI] += diag()[faceLabels[edgeI]];

    boundaryCoeffs_[patchI][edgeI] = value;
}


template<>
Foam::solverPerformance Foam::faMatrix<Foam::scalar>::solve
(
    const dictionary& solverControls
)
{
    DebugInFunction
        << "solving faMatrix<scalar>"
        << endl;

    const int logLevel =
        solverControls.getOrDefault<int>
        (
            "log",
            solverPerformance::debug
        );

    auto& psi =
        const_cast<GeometricField<scalar, faPatchField, areaMesh>&>(psi_);

    scalarField saveDiag(diag());
    addBoundaryDiag(diag(), 0);

    scalarField totalSource(source_);
    addBoundarySource(totalSource, 0);

    // Solver call
    solverPerformance solverPerf = lduMatrix::solver::New
    (
        psi_.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        psi_.boundaryField().scalarInterfaces(),
        solverControls
    )->solve(psi.ref(), totalSource);

    if (logLevel)
    {
        solverPerf.print(Info);
    }

    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi.mesh().data().setSolverPerformance(psi.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::tmp<Foam::scalarField> Foam::faMatrix<Foam::scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), Zero);
    addBoundaryDiag(boundaryDiag, 0);

    const scalarField& psif = psi_.internalField();
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
Foam::tmp<Foam::areaScalarField> Foam::faMatrix<Foam::scalar>::H() const
{
    auto tHphi = areaScalarField::New
    (
        "H(" + psi_.name() + ')',
        psi_.mesh(),
        dimensions_/dimArea,
        faPatchFieldBase::extrapolatedCalculatedType()
    );
    auto& Hphi = tHphi.ref();

    Hphi.primitiveFieldRef() = (lduMatrix::H(psi_.primitiveField()) + source_);
    addBoundarySource(Hphi.primitiveFieldRef());

    Hphi.ref() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// ************************************************************************* //
