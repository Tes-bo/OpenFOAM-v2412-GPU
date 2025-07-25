/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2016 OpenFOAM Foundation
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

#include "noLift.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace liftModels
{
    defineTypeNameAndDebug(noLift, 0);
    addToRunTimeSelectionTable(liftModel, noLift, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::liftModels::noLift::noLift
(
    const dictionary& dict,
    const phasePair& pair
)
:
    liftModel(dict, pair)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::liftModels::noLift::~noLift()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::liftModels::noLift::Cl() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();

    return volScalarField::New
    (
        "Cl",
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimless, Zero)
    );
}


Foam::tmp<Foam::volVectorField> Foam::liftModels::noLift::F() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();

    return volVectorField::New
    (
        IOobject::scopedName("noLift", "F"),
        IOobject::NO_REGISTER,
        mesh,
        dimensionedVector(dimF, Zero)
    );
}


Foam::tmp<Foam::surfaceScalarField> Foam::liftModels::noLift::Ff() const
{
    const fvMesh& mesh = this->pair_.phase1().mesh();

    return surfaceScalarField::New
    (
        IOobject::scopedName("noLift", "Ff"),
        IOobject::NO_REGISTER,
        mesh,
        dimensionedScalar(dimF*dimArea, Zero)
    );
}


// ************************************************************************* //
