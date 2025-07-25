/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2012 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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

#include "solidThermo.H"
#include "fvMesh.H"


/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(solidThermo, 0);
    defineRunTimeSelectionTable(solidThermo, fvMesh);
    defineRunTimeSelectionTable(solidThermo, dictionary);
    defineRunTimeSelectionTable(solidThermo, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    rhoThermo(mesh, phaseName)
{}


Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    rhoThermo(mesh, dict, phaseName)
{}


Foam::solidThermo::solidThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
     rhoThermo(mesh, phaseName, dictionaryName)
{}



// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(mesh, phaseName);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
{
    return basicThermo::New<solidThermo>(mesh, dict, phaseName);
}


Foam::autoPtr<Foam::solidThermo> Foam::solidThermo::New
(
     const fvMesh& mesh,
     const word& phaseName,
     const word& dictName
)
{
    return basicThermo::New<solidThermo>(mesh, phaseName, dictName);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::solidThermo::~solidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField> Foam::solidThermo::rho() const
{
    return rho_;
}


Foam::tmp<Foam::scalarField> Foam::solidThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}


Foam::volScalarField& Foam::solidThermo::rho()
{
    return rho_;
}


bool Foam::solidThermo::read()
{
    return regIOobject::read();
}


// ************************************************************************* //
