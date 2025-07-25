/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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

#include "noBlending.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace blendingMethods
{
    defineTypeNameAndDebug(noBlending, 0);

    addToRunTimeSelectionTable
    (
        blendingMethod,
        noBlending,
        dictionary
    );
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::noBlending
(
    const dictionary& dict,
    const wordList& phaseNames
)
:
    blendingMethod(dict),
    continuousPhase_(dict.get<word>("continuousPhase"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::blendingMethods::noBlending::~noBlending()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f1
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh = phase1.mesh();

    return volScalarField::New
    (
        "f",
        IOobject::NO_REGISTER,
        mesh,
        scalar(phase2.name() != continuousPhase_),  // value 0/1
        dimless
    );
}


Foam::tmp<Foam::volScalarField> Foam::blendingMethods::noBlending::f2
(
    const phaseModel& phase1,
    const phaseModel& phase2
) const
{
    const fvMesh& mesh = phase1.mesh();

    return volScalarField::New
    (
        "f",
        IOobject::NO_REGISTER,
        mesh,
        scalar(phase1.name() == continuousPhase_),  // value 0/1
        dimless
    );
}


// ************************************************************************* //
