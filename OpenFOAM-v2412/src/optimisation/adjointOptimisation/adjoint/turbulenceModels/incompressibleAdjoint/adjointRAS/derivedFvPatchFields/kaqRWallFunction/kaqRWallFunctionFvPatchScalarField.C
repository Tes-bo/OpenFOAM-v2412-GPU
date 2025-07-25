/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2022 PCOpt/NTUA
    Copyright (C) 2014-2022 FOSS GP
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

#include "kaqRWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "nutkWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::kaqRWallFunctionFvPatchScalarField::kaqRWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<scalar>(p, iF),
    adjointScalarBoundaryCondition(p, iF, word::null)
{}


Foam::kaqRWallFunctionFvPatchScalarField::kaqRWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    kqRWallFunctionFvPatchField<scalar>(p, iF, dict),
    adjointScalarBoundaryCondition(p, iF, dict.get<word>("solverName"))
{}


Foam::kaqRWallFunctionFvPatchScalarField::kaqRWallFunctionFvPatchScalarField
(
    const kaqRWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    kqRWallFunctionFvPatchField<scalar>(ptf, p, iF, mapper),
    adjointScalarBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::kaqRWallFunctionFvPatchScalarField::kaqRWallFunctionFvPatchScalarField
(
    const kaqRWallFunctionFvPatchScalarField& tkqrwfpf
)
:
    kqRWallFunctionFvPatchField<scalar>(tkqrwfpf),
    adjointScalarBoundaryCondition(tkqrwfpf)
{}


Foam::kaqRWallFunctionFvPatchScalarField::kaqRWallFunctionFvPatchScalarField
(
    const kaqRWallFunctionFvPatchScalarField& tkqrwfpf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    kqRWallFunctionFvPatchField<scalar>(tkqrwfpf, iF),
    adjointScalarBoundaryCondition(tkqrwfpf)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::kaqRWallFunctionFvPatchScalarField::manipulateMatrix
(
    fvMatrix<scalar>& matrix
)
{
    scalarField& source = matrix.source();
    tmp<fvPatchScalarField> tnutWall(boundaryContrPtr_->turbulentDiffusivity());
    fvPatchScalarField& nutWall = tnutWall.constCast();
    if (isA<nutkWallFunctionFvPatchScalarField>(nutWall))
    {
        const label patchi(patch().index());
        const scalarField& magSf = patch().magSf();
        const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
        (
            IOobject::groupName
            (
                turbulenceModel::propertiesName,
                internalField().group()
            )
        );
        const scalarField& y = turbModel.y()[patchi];
        const tmp<scalarField> tnuw = turbModel.nu(patchi);
        const scalarField& nuw = tnuw();
        const nutWallFunctionFvPatchScalarField& nutWF =
            refCast<nutWallFunctionFvPatchScalarField>(nutWall);
        const wallFunctionCoefficients& wallCoeffs = nutWF.wallCoeffs();
        const scalar Cmu = wallCoeffs.Cmu();
        const scalar kappa = wallCoeffs.kappa();
        const scalar E = wallCoeffs.E();
        const scalar yPlusLam = wallCoeffs.yPlusLam();

        const scalar Cmu25 = pow025(Cmu);

        const labelList& faceCells = patch().faceCells();
        const fvPatchVectorField& Uw = boundaryContrPtr_->Ub();
        const scalarField magGradUw(mag(Uw.snGrad()));

        tmp<scalarField> tdJdnut = boundaryContrPtr_->dJdnut();
        const scalarField& dJdnut = tdJdnut();

        const tmp<volScalarField> tk = turbModel.k();
        const volScalarField& k = tk();
        forAll(dJdnut, facei)
        {
            const label celli = faceCells[facei];

            const scalar sqrtkCell(sqrt(k[celli]));
            const scalar yPlus = Cmu25*y[facei]*sqrtkCell/nuw[facei];
            if (yPlusLam < yPlus)
            {
                const scalar logEyPlus = log(E*yPlus);
                const scalar dnut_dyPlus =
                        nuw[facei]*kappa*(logEyPlus - 1)/sqr(logEyPlus);
                const scalar dyPlus_dk =
                    Cmu25*y[facei]/(2*nuw[facei]*sqrtkCell);
                const scalar dnut_dk = dnut_dyPlus*dyPlus_dk;
                source[celli] -= dJdnut[facei]*dnut_dk*magSf[facei];
            }
        }
    }
}


void Foam::kaqRWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    kqRWallFunctionFvPatchField<scalar>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        kaqRWallFunctionFvPatchScalarField
    );
}


// ************************************************************************* //
