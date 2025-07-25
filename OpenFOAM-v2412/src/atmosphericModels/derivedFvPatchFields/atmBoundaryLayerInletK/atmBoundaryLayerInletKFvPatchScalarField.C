/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "atmBoundaryLayerInletKFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayer(iF.time(), p.patch())
{}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    atmBoundaryLayer(iF.time(), p.patch(), dict)
{
    phiName_ = dict.getOrDefault<word>("phi", "phi");

    refValue() = k(patch().Cf());
    refGrad() = 0;
    valueFraction() = 1;

    if (!initABL_)
    {
        this->readValueEntry(dict, IOobjectOption::MUST_READ);
    }
    else
    {
        scalarField::operator=(refValue());
        initABL_ = false;
    }
}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const atmBoundaryLayerInletKFvPatchScalarField& psf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(psf, p, iF, mapper),
    atmBoundaryLayer(psf, p, mapper)
{}


atmBoundaryLayerInletKFvPatchScalarField::
atmBoundaryLayerInletKFvPatchScalarField
(
    const atmBoundaryLayerInletKFvPatchScalarField& psf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(psf, iF),
    atmBoundaryLayer(psf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void atmBoundaryLayerInletKFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    refValue() = k(patch().Cf());

    inletOutletFvPatchScalarField::updateCoeffs();
}


void atmBoundaryLayerInletKFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);
    atmBoundaryLayer::autoMap(m);
}


void atmBoundaryLayerInletKFvPatchScalarField::rmap
(
    const fvPatchScalarField& psf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(psf, addr);

    const atmBoundaryLayerInletKFvPatchScalarField& blpsf =
        refCast<const atmBoundaryLayerInletKFvPatchScalarField>(psf);

    atmBoundaryLayer::rmap(blpsf, addr);
}


void atmBoundaryLayerInletKFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    os.writeEntryIfDifferent<word>("phi", "phi", phiName_);
    atmBoundaryLayer::write(os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    atmBoundaryLayerInletKFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
