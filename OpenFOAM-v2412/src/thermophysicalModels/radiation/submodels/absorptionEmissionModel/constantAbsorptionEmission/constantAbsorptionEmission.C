/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "constantAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(constantAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            constantAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::constantAbsorptionEmission::constantAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
    a_("absorptivity", coeffsDict_),
    e_("emissivity", coeffsDict_),
    E_("E", coeffsDict_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::constantAbsorptionEmission::aCont(const label bandI) const
{
    return volScalarField::New
    (
        "a",
        IOobject::NO_REGISTER,
        mesh_,
        a_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::constantAbsorptionEmission::eCont(const label bandI) const
{
    return volScalarField::New
    (
        "e",
        IOobject::NO_REGISTER,
        mesh_,
        e_
    );
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::constantAbsorptionEmission::ECont(const label bandI) const
{
    return volScalarField::New
    (
        "E",
        IOobject::NO_REGISTER,
        mesh_,
        E_
    );
}


// ************************************************************************* //
