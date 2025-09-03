/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd
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

#include "alphatJayatillekeWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "nutWallFunctionFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

scalar alphatJayatillekeWallFunctionFvPatchScalarField::tolerance_ = 0.01;
label alphatJayatillekeWallFunctionFvPatchScalarField::maxIters_ = 10;

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::checkType()
{
    if (!isA<wallFvPatch>(patch()))
    {
        FatalErrorInFunction
            << "Invalid wall function specification" << nl
            << "    Patch type for patch " << patch().name()
            << " must be wall" << nl
            << "    Current patch type is " << patch().type() << nl << endl
            << abort(FatalError);
    }
}


tmp<scalarField> alphatJayatillekeWallFunctionFvPatchScalarField::yPlus
(
    const turbulenceModel& turbModel
) const
{
    const label patchi = patch().index();
    const tmp<volScalarField> tnut = turbModel.nut();
    const volScalarField& nut = tnut();

    if (isA<nutWallFunctionFvPatchScalarField>(nut.boundaryField()[patchi]))
    {
        const auto& nutPf =
            dynamic_cast<const nutWallFunctionFvPatchScalarField&>
            (
                nut.boundaryField()[patchi]
            );

        return nutPf.yPlus();
    }
    else
    {
        const scalarField& y = turbModel.y()[patchi];
        const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];

        return
            y*sqrt(turbModel.nuEff(patchi)*mag(Uw.snGrad()))
           /turbModel.nu(patchi);
    }
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::Psmooth
(
    const scalar Prat
) const
{
    return 9.24*(pow(Prat, 0.75) - 1.0)*(1.0 + 0.28*exp(-0.007*Prat));
}


scalar alphatJayatillekeWallFunctionFvPatchScalarField::yPlusTherm
(
    const scalar P,
    const scalar Prat
) const
{
    scalar ypt = 11;

    for (int iter = 0; iter < maxIters_; ++iter)
    {
        const scalar f = ypt - (log(E_*ypt)/kappa_ + P)/Prat;
        const scalar df = 1.0 - 1.0/(ypt*kappa_*Prat);
        const scalar yptNew = ypt - f/df;

        if (yptNew < VSMALL)
        {
            return 0;
        }
        else if (mag(yptNew - ypt) < tolerance_)
        {
            return yptNew;
        }
        else
        {
            ypt = yptNew;
        }
     }

    return ypt;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Prt_(0.85),
    kappa_(0.41),
    E_(9.8)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Prt_(ptf.Prt_),
    kappa_(ptf.kappa_),
    E_(ptf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    Prt_(dict.get<scalar>("Prt")), // force read to avoid ambiguity
    kappa_(dict.getOrDefault<scalar>("kappa", 0.41)),
    E_(dict.getOrDefault<scalar>("E", 9.8))
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    Prt_(wfpsf.Prt_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{
    checkType();
}


alphatJayatillekeWallFunctionFvPatchScalarField::
alphatJayatillekeWallFunctionFvPatchScalarField
(
    const alphatJayatillekeWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    Prt_(wfpsf.Prt_),
    kappa_(wfpsf.kappa_),
    E_(wfpsf.E_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatJayatillekeWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    // Retrieve turbulence properties from model

    const auto& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const scalarField yPlusp(yPlus(turbModel));

    const tmp<volScalarField> tnu = turbModel.nu();
    const volScalarField& nu = tnu();
    const scalarField& nuw = nu.boundaryField()[patchi];

    const auto& transportProperties =
        db().lookupObject<IOdictionary>("transportProperties");

    // Molecular Prandtl number
    const scalar Pr
    (
        dimensionedScalar("Pr", dimless, transportProperties).value()
    );

    // Populate boundary values
    scalarField& alphatw = *this;
    forAll(alphatw, facei)
    {
        const scalar yPlus = yPlusp[facei];

        // Molecular-to-turbulent Prandtl number ratio
        const scalar Prat = Pr/Prt_;

        // Thermal sublayer thickness
        const scalar P = Psmooth(Prat);
        const scalar yPlusTherm = this->yPlusTherm(P, Prat);

        // Update turbulent thermal conductivity
        if (yPlus > yPlusTherm)
        {
            const scalar nu = nuw[facei];
            const scalar kt =
                nu*(yPlus/(Prt_*(log(E_*yPlus)/kappa_ + P)) - 1.0/Pr);

            alphatw[facei] = max(scalar(0), kt);
        }
        else
        {
            alphatw[facei] = 0.0;
        }
    }

    fixedValueFvPatchField<scalar>::updateCoeffs();
}


void alphatJayatillekeWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntry("Prt", Prt_);
    os.writeEntryIfDifferent<scalar>("kappa", 0.41, kappa_);
    os.writeEntryIfDifferent<scalar>("E", 9.8, E_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    alphatJayatillekeWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //
