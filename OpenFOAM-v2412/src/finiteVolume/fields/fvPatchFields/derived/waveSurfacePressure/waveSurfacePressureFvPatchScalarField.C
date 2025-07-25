/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "waveSurfacePressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "gravityMeshObject.H"
#include "EulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "backwardDdtScheme.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::waveSurfacePressureFvPatchScalarField::ddtSchemeType
>
Foam::waveSurfacePressureFvPatchScalarField::ddtSchemeTypeNames_
({
    {
        ddtSchemeType::tsEuler,
        fv::EulerDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsCrankNicolson,
        fv::CrankNicolsonDdtScheme<scalar>::typeName_()
    },
    {
        ddtSchemeType::tsBackward,
        fv::backwardDdtScheme<scalar>::typeName_()
    },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    phiName_("phi"),
    zetaName_("zeta"),
    rhoName_("rho")
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    zetaName_(dict.getOrDefault<word>("zeta", "zeta")),
    rhoName_(dict.getOrDefault<word>("rho", "rho"))
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    phiName_(ptf.phiName_),
    zetaName_(ptf.zetaName_),
    rhoName_(ptf.rhoName_)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf
)
:
    fixedValueFvPatchScalarField(wspsf),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_)
{}


Foam::waveSurfacePressureFvPatchScalarField::
waveSurfacePressureFvPatchScalarField
(
    const waveSurfacePressureFvPatchScalarField& wspsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wspsf, iF),
    phiName_(wspsf.phiName_),
    zetaName_(wspsf.zetaName_),
    rhoName_(wspsf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::waveSurfacePressureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const label patchi = patch().index();

    const scalar dt = db().time().deltaTValue();

    // Retrieve non-const access to zeta field from the database
    volVectorField& zeta = db().lookupObjectRef<volVectorField>(zetaName_);
    vectorField& zetap = zeta.boundaryFieldRef()[patchi];

    // Lookup d/dt scheme from database for zeta
    const word ddtSchemeName(zeta.mesh().ddtScheme(zeta.name()));
    ddtSchemeType ddtScheme(ddtSchemeTypeNames_[ddtSchemeName]);

    // Retrieve the flux field from the database
    const auto& phi = db().lookupObject<surfaceScalarField>(phiName_);

    // Cache the patch face-normal vectors
    tmp<vectorField> nf(patch().nf());

    // Change in zeta due to flux
    vectorField dZetap(dt*nf()*phi.boundaryField()[patchi]/patch().magSf());

    if (phi.dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);

        dZetap /= rhop;
    }

    const volVectorField& zeta0 = zeta.oldTime();

    switch (ddtScheme)
    {
        case tsEuler:
        case tsCrankNicolson:
        {
            zetap = zeta0.boundaryField()[patchi] + dZetap;

            break;
        }
        case tsBackward:
        {
            scalar dt0 = db().time().deltaT0Value();

            scalar c = 1.0 + dt/(dt + dt0);
            scalar c00 = dt*dt/(dt0*(dt + dt0));
            scalar c0 = c + c00;

            zetap =
                (
                    c0*zeta0.boundaryField()[patchi]
                  - c00*zeta0.oldTime().boundaryField()[patchi]
                  + dZetap
                )/c;

            break;
        }
        default:
        {
            FatalErrorInFunction
                << ddtSchemeName << nl
                << "    on patch " << this->patch().name()
                << " of field " << this->internalField().name()
                << " in file " << this->internalField().objectPath()
                << abort(FatalError);
        }
    }


    Info<< "min/max zetap = " << gMin(zetap & nf()) << ", "
        << gMax(zetap & nf()) << endl;

    // Update the surface pressure
    const uniformDimensionedVectorField& g =
        meshObjects::gravity::New(db().time());

    operator==(-g.value() & zetap);

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::waveSurfacePressureFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("zeta", "zeta", zetaName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        waveSurfacePressureFvPatchScalarField
    );
}

// ************************************************************************* //
