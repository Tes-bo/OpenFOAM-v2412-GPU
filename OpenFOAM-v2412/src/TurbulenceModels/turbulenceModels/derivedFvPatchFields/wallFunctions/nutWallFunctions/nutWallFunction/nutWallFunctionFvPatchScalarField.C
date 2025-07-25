/*---------------------------------------------------------------------------* \
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016, 2019 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "nutWallFunctionFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "turbulenceModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(nutWallFunctionFvPatchScalarField, 0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::nutWallFunctionFvPatchScalarField::checkType()
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


const Foam::volVectorField& Foam::nutWallFunctionFvPatchScalarField::U
(
    const turbulenceModel& turb
) const
{
    if (UName_.empty())
    {
        return turb.U();
    }

    return db().lookupObject<volVectorField>(UName_);
}


void Foam::nutWallFunctionFvPatchScalarField::writeLocalEntries
(
    Ostream& os
) const
{
    os.writeEntryIfDifferent<word>("U", word::null, UName_);
    wallCoeffs_.writeEntries(os);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    UName_(),
    wallCoeffs_()
{
    checkType();
}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    UName_(ptf.UName_),
    wallCoeffs_(ptf.wallCoeffs_)
{
    checkType();
}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    UName_(dict.getOrDefault<word>("U", word::null)),
    wallCoeffs_(dict)
{
    checkType();
}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& wfpsf
)
:
    fixedValueFvPatchScalarField(wfpsf),
    UName_(wfpsf.UName_),
    wallCoeffs_(wfpsf.wallCoeffs_)
{
    checkType();
}


Foam::nutWallFunctionFvPatchScalarField::nutWallFunctionFvPatchScalarField
(
    const nutWallFunctionFvPatchScalarField& wfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(wfpsf, iF),
    UName_(wfpsf.UName_),
    wallCoeffs_(wfpsf.wallCoeffs_)
{
    checkType();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::nutWallFunctionFvPatchScalarField&
Foam::nutWallFunctionFvPatchScalarField::nutw
(
    const turbulenceModel& turbModel,
    const label patchi
)
{
    return
        refCast<const nutWallFunctionFvPatchScalarField>
        (
            turbModel.nut()().boundaryField()[patchi],
            patchi
        );
}


void Foam::nutWallFunctionFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    operator==(calcNut());

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::nutWallFunctionFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
}


// ************************************************************************* //
