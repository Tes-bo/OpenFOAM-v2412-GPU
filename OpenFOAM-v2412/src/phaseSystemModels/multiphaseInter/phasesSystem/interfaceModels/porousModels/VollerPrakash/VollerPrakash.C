/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2020 OpenCFD Ltd.
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

#include "VollerPrakash.H"
#include "phasePair.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace porousModels
{
    defineTypeNameAndDebug(VollerPrakash, 0);
    addToRunTimeSelectionTable
    (
        porousModel,
        VollerPrakash,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::porousModels::VollerPrakash::VollerPrakash
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    porousModel(dict, mesh),
    Cu_(dict.get<scalar>("Cu")),
    solidPhase_(dict.get<word>("solidPhase"))
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::porousModels::VollerPrakash::S() const
{
    const volScalarField& solidAlpha =
        mesh_.lookupObject<volScalarField>(solidPhase_);

    if (debug && mesh_.time().writeTime())
    {
        volScalarField Su
        (
            "Su" + solidAlpha.name(),
            Cu_*sqr(solidAlpha)/(pow3(1.0 - solidAlpha) + 1e-3)
        );
        Su.write();
    }

    return Cu_*sqr(solidAlpha)/(pow3(1.0 - solidAlpha) + 1e-3);
}


// ************************************************************************* //
