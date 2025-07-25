/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "adjointInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointInletVelocityFvPatchVectorField::
adjointInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, word::null)
{}


Foam::adjointInletVelocityFvPatchVectorField::
adjointInletVelocityFvPatchVectorField
(
    const adjointInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    adjointVectorBoundaryCondition(p, iF, ptf.adjointSolverName_)
{}


Foam::adjointInletVelocityFvPatchVectorField::
adjointInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    adjointVectorBoundaryCondition(p, iF, dict.get<word>("solverName"))
{
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::adjointInletVelocityFvPatchVectorField::
adjointInletVelocityFvPatchVectorField
(
    const adjointInletVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    adjointVectorBoundaryCondition(pivpvf)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::adjointInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    // Objective function contribution
    tmp<vectorField> tsource = boundaryContrPtr_->normalVelocitySource();
    const vectorField& source = tsource();

    operator==(-source);

    fixedValueFvPatchVectorField::updateCoeffs();
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointInletVelocityFvPatchVectorField::valueInternalCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector>>::New(this->size(), pTraits<vector>::one);
}


Foam::tmp<Foam::Field<Foam::vector>>
Foam::adjointInletVelocityFvPatchVectorField::valueBoundaryCoeffs
(
    const tmp<scalarField>&
) const
{
    return tmp<Field<vector>>::New(this->size(), Zero);
}


void Foam::adjointInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchField<vector>::write(os);
    os.writeEntry("solverName", adjointSolverName_);
    fvPatchField<vector>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
