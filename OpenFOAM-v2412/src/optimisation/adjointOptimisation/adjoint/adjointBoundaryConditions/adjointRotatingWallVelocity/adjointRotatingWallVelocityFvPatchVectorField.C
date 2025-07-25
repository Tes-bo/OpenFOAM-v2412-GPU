/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 PCOpt/NTUA
    Copyright (C) 2020 FOSS GP
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

#include "adjointRotatingWallVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::adjointRotatingWallVelocityFvPatchVectorField::
adjointRotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    adjointWallVelocityFvPatchVectorField(p, iF),
    origin_(),
    axis_(Zero),
    omega_(nullptr)
{}


Foam::adjointRotatingWallVelocityFvPatchVectorField::
adjointRotatingWallVelocityFvPatchVectorField
(
    const adjointRotatingWallVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    adjointWallVelocityFvPatchVectorField(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_.clone())
{}


Foam::adjointRotatingWallVelocityFvPatchVectorField::
adjointRotatingWallVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    adjointWallVelocityFvPatchVectorField(p, iF, dict),
    origin_(dict.get<vector>("origin")),
    axis_(dict.get<vector>("axis")),
    omega_(Function1<scalar>::New("omega", dict, &db()))
{}


Foam::adjointRotatingWallVelocityFvPatchVectorField::
adjointRotatingWallVelocityFvPatchVectorField
(
    const adjointRotatingWallVelocityFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    adjointWallVelocityFvPatchVectorField(pivpvf, iF),
    origin_(pivpvf.origin_),
    axis_(pivpvf.axis_),
    omega_(pivpvf.omega_.clone())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::tensorField>
Foam::adjointRotatingWallVelocityFvPatchVectorField::dxdbMult() const
{
    const scalar t(this->db().time().timeOutputValue());
    const scalar om(omega_->value(t));
    const vector omega(om*axis_/mag(axis_));
    tensor mult
    (
        scalar(0), -omega.z(),  omega.y(),
        omega.z(),  scalar(0), -omega.x(),
       -omega.y(),  omega.x(),  scalar(0)
    );

    return tmp<tensorField>::New(patch().size(), mult);
}


void Foam::adjointRotatingWallVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    adjointWallVelocityFvPatchVectorField::write(os);
    os.writeEntry("origin", origin_);
    os.writeEntry("axis", axis_);
    omega_->writeData(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        adjointRotatingWallVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
