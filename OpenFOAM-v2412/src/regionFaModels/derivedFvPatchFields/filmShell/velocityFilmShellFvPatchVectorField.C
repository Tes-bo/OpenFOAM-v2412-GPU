/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "velocityFilmShellFvPatchVectorField.H"
#include "dictionaryContent.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(p, iF),
    baffle_(nullptr),
    dict_(),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = 1;
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchField<vector>
    (
        ptf,
        p,
        iF,
        mapper
    ),
    baffle_(nullptr),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchField<vector>(p, iF),
    baffle_(nullptr),
    dict_
    (
        // Copy dictionary, but without "heavy" data chunks
        dictionaryContent::copyDict
        (
            dict,
            wordList(),  // allow
            wordList     // deny
            ({
                "type",  // redundant
                "value", "refValue", "refGradient", "valueFraction"
            })
        )
    ),
    curTimeIndex_(-1),
    zeroWallVelocity_(dict.getOrDefault<bool>("zeroWallVelocity", true))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);

    if (this->readMixedEntries(dict))
    {
        // Full restart
    }
    else
    {
        // Start from user entered data. Assume fixedValue.
        refValue() = *this;
        refGrad() = Zero;
        valueFraction() = 1;
    }

    if (!baffle_)
    {
        baffle_.reset(baffleType::New(p.boundaryMesh().mesh(), dict_));
    }
}


velocityFilmShellFvPatchVectorField::velocityFilmShellFvPatchVectorField
(
    const velocityFilmShellFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    mixedFvPatchField<vector>(ptf, iF),
    baffle_(nullptr),
    dict_(ptf.dict_),
    curTimeIndex_(-1),
    zeroWallVelocity_(true)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void velocityFilmShellFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Execute the change only once per time-step
    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        baffle_->evolve();

        vectorField& pfld = *this;

        baffle_->vsm().mapToVolumePatch(baffle_->Us(), pfld, patch().index());

        refGrad() = Zero;
        valueFraction() = 1;

        if (zeroWallVelocity_)
        {
            refValue() = Zero;
        }
        else
        {
            refValue() = pfld;
        }

        curTimeIndex_ = this->db().time().timeIndex();
    }

    mixedFvPatchField<vector>::updateCoeffs();
}


void velocityFilmShellFvPatchVectorField::write(Ostream& os) const
{
    mixedFvPatchField<vector>::write(os);
    dict_.write(os, false);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    velocityFilmShellFvPatchVectorField
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
