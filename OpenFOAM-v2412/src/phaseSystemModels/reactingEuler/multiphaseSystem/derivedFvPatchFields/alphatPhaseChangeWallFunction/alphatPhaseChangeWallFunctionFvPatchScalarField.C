/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2018 OpenFOAM Foundation
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

#include "alphatPhaseChangeWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(alphatPhaseChangeWallFunctionFvPatchScalarField, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    dmdt_(p.size(), Zero),
    mDotL_(p.size(), Zero)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    dmdt_("dmdt", dict, p.size(), IOobjectOption::LAZY_READ),
    mDotL_("mDotL", dict, p.size(), IOobjectOption::LAZY_READ)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    dmdt_(ptf.dmdt_, mapper),
    mDotL_(ptf.mDotL_, mapper)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& awfpsf
)
:
    fixedValueFvPatchScalarField(awfpsf),
    dmdt_(awfpsf.dmdt_),
    mDotL_(awfpsf.mDotL_)
{}


alphatPhaseChangeWallFunctionFvPatchScalarField::
alphatPhaseChangeWallFunctionFvPatchScalarField
(
    const alphatPhaseChangeWallFunctionFvPatchScalarField& awfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(awfpsf, iF),
    dmdt_(awfpsf.dmdt_),
    mDotL_(awfpsf.mDotL_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void alphatPhaseChangeWallFunctionFvPatchScalarField::
write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    dmdt_.writeEntry("dmdt", os);
    mDotL_.writeEntry("mDotL", os);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible
} // End namespace Foam

// ************************************************************************* //
