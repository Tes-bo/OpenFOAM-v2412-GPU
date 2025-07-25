/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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

#include "adjointWallVelocityFvPatchVectorField.H"
#include "nutUSpaldingWallFunctionFvPatchScalarField.H"
#include "fvMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointWallVelocityFvPatchVectorField::
adjointWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, word::null),
    kappa_(0.41),
    E_(9.8)
{}


Foam::adjointWallVelocityFvPatchVectorField::
adjointWallVelocityFvPatchVectorField
(
    const adjointWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointVectorBoundaryCondition(p, iF, ptf.adjointSolverName_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{}


Foam::adjointWallVelocityFvPatchVectorField::
adjointWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, dict.get<word>("solverName")),
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    E_(dict.getOrDefault<scalar>("E", 9.8))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::adjointWallVelocityFvPatchVectorField::
adjointWallVelocityFvPatchVectorField
(
    const adjointWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointVectorBoundaryCondition(pivpvf),
    kappa_(pivpvf.kappa_),
    E_(pivpvf.E_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointWallVelocityFvPatchVectorField::manipulateMatrix
(
    fvMatrix<vector>& matrix
)
{
    addProfiling
    (
        adjointWallVelocityFvPatchVectorField,
        "adjointWallVelocityFvPatchVectorField::manipulateMatrix"
    );
    // Grab ref to the diagonal matrix
    vectorField& source = matrix.source();

    // Define boundary condition type
    typedef Foam::nutUSpaldingWallFunctionFvPatchScalarField
        SAwallFunctionPatchField;

    tmp<fvPatchScalarField> nutPatch(boundaryContrPtr_->turbulentDiffusivity());
    if
    (
        isA<SAwallFunctionPatchField>(nutPatch())
     && patch().size() != 0
    )
    {
        const tmp<vectorField> tnf = patch().nf();
        const vectorField& nf = tnf();
        const scalarField& magSf = patch().magSf();

        const fvPatchField<vector>& Up = boundaryContrPtr_->Ub();
        const fvPatchField<vector>& Uap = *this;

        const vectorField Uc(Up.patchInternalField());
        const vectorField Uc_t(Uc - (Uc & nf)*nf);

        // By convention, tf has the direction of the tangent PRIMAL velocity
        // at the first cell off the wall
        const vectorField tf(Uc_t/mag(Uc_t));

        tmp<scalarField> tnuw = boundaryContrPtr_->momentumDiffusion();
        const scalarField& nuw = tnuw();
        tmp<scalarField> tnu = boundaryContrPtr_->laminarDiffusivity();
        const scalarField& nu = tnu();
        tmp<scalarField> tyC = boundaryContrPtr_->wallDistance();
        const scalarField& yC = tyC();

        const scalarField magGradU(mag(Up.snGrad()));
        const scalarField vtau(sqrt(nuw*magGradU));
        const scalarField uPlus(mag(Uc)/vtau);
        const scalarField yPlus(yC*vtau/nu);
        const scalarField kUu(min(kappa_*uPlus, scalar(50)));
        const scalarField auxA((kappa_/E_)*(exp(kUu)-1 - kUu - 0.5*kUu*kUu));
        const scalarField auxB(-(1 + auxA)/(yPlus + uPlus*(1 + auxA)));

        // Tangent components are according to tf
        tmp<vectorField> tsource = boundaryContrPtr_->normalVelocitySource();
        const scalarField rt(tsource() & tf);
        const scalarField Uap_t(Uap & tf);

        forAll(Up, faceI)
        {
            label cellI = patch().faceCells()[faceI];
            source[cellI] +=
                2*auxB[faceI]*vtau[faceI]*((rt[faceI] + Uap_t[faceI]))
               *(Uc[faceI]/mag(Uc[faceI]))*magSf[faceI];
        }
    }
}


void Foam::adjointWallVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const fvPatchField<vector>& Up = boundaryContrPtr_->Ub();

    // Patch geometry
    tmp<vectorField> tnf = patch().nf();
    const vectorField& nf = tnf();

    // Internal fields
    vectorField Uac(this->patchInternalField());
    vectorField Uc(Up.patchInternalField());

    // Tangent vector based on the  direction of Vc
    vectorField Uc_t(Uc - (Uc & nf)*nf);
    vectorField tf1(Uc_t/mag(Uc_t));

    // Tangent vector as the cross product of tf1 x nf
    vectorField tf2((tf1 ^ nf)/mag(tf1 ^ nf));

    // Normal adjoint component comes from the objective function
    tmp<vectorField> tsource = boundaryContrPtr_->normalVelocitySource();
    vectorField Uan(-(tsource() & nf)*nf);

    // Tangential adjoint velocity in the t1 direction depends on the primal
    // wall function used
    vectorField Uap_t1(patch().size(), Zero);
    typedef Foam::nutUSpaldingWallFunctionFvPatchScalarField
        SAwallFunctionPatchField;

    tmp<fvPatchScalarField> nutb(boundaryContrPtr_->turbulentDiffusivity());
    if (isA<SAwallFunctionPatchField>(nutb()))
    {
        Uap_t1 = (Uac & tf1)*tf1;
        // leaving out second term for now
        //- (1./delta)*((gradUaC & nf) & tf1)*tf1;
    }
    else
    {
        Uap_t1 = - (tsource() & tf1)*tf1;
    }

    // Adjoint velocity in the t2 direction
    vectorField Uap_t2(-(tsource() & tf2)*tf2);

    operator==(Uan + Uap_t1 + Uap_t2);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::adjointWallVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("kappa", kappa_);
    os.writeEntry("E", E_);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
