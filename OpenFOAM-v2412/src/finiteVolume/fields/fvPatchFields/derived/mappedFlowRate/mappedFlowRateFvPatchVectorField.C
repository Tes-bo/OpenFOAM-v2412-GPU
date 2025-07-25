/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "mappedFlowRateFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "mappedPatchBase.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedFlowRateFvPatchVectorField::mappedFlowRateFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    nbrPhiName_("none"),
    phiName_("phi"),
    rhoName_("rho")
{}


Foam::mappedFlowRateFvPatchVectorField::mappedFlowRateFvPatchVectorField
(
    const mappedFlowRateFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::mappedFlowRateFvPatchVectorField::mappedFlowRateFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict),
    nbrPhiName_(dict.getOrDefault<word>("nbrPhi", "phi")),
    phiName_(dict.getOrDefault<word>("phi", "phi")),
    rhoName_(dict.getOrDefault<word>("rho", "rho"))
{}


Foam::mappedFlowRateFvPatchVectorField::mappedFlowRateFvPatchVectorField
(
    const mappedFlowRateFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


Foam::mappedFlowRateFvPatchVectorField::mappedFlowRateFvPatchVectorField
(
    const mappedFlowRateFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF),
    nbrPhiName_(ptf.nbrPhiName_),
    phiName_(ptf.phiName_),
    rhoName_(ptf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mappedFlowRateFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Since we're inside initEvaluate/evaluate there might be processor
    // comms underway. Change the tag we use.
    const int oldTag = UPstream::incrMsgType();

    // Get the coupling information from the mappedPatchBase
    const mappedPatchBase& mpp = refCast<const mappedPatchBase>
    (
        patch().patch()
    );
    const polyMesh& nbrMesh = mpp.sampleMesh();
    const fvPatch& nbrPatch = refCast<const fvMesh>
    (
        nbrMesh
    ).boundary()[mpp.samplePolyPatch().index()];

    scalarList phi
    (
        nbrPatch.lookupPatchField<surfaceScalarField>(nbrPhiName_)
    );

    mpp.distribute(phi);


    const surfaceScalarField& phiName =
        db().lookupObject<surfaceScalarField>(phiName_);

    scalarField U(-phi/patch().magSf());

    vectorField n(patch().nf());

    if (phiName.dimensions() == dimVolume/dimTime)
    {
        // volumetric flow-rate
        operator==(n*U);
    }
    else if (phiName.dimensions() == dimMass/dimTime)
    {
        const auto& rhop = patch().lookupPatchField<volScalarField>(rhoName_);

        // mass flow-rate
        operator==(n*U/rhop);

        if (debug)
        {
            scalar phi = gSum(rhop*(*this) & patch().Sf());
            Info<< patch().boundaryMesh().mesh().name() << ':'
                << patch().name() << ':'
                << this->internalField().name() << " <- "
                << nbrMesh.name() << ':'
                << nbrPatch.name() << ':'
                << this->internalField().name() << " :"
                << " mass flux[Kg/s]:" << -phi
                << endl;
        }
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of " << phiName_ << " are incorrect" << nl
            << "    on patch " << this->patch().name()
            << " of field " << this->internalField().name()
            << " in file " << this->internalField().objectPath()
            << nl << exit(FatalError);
    }

    UPstream::msgType(oldTag);  // Restore tag

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::mappedFlowRateFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchField<vector>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    os.writeEntryIfDifferent<word>("rho", "rho", rhoName_);
    os.writeEntry("nbrPhi", nbrPhiName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       mappedFlowRateFvPatchVectorField
   );
}


// ************************************************************************* //
