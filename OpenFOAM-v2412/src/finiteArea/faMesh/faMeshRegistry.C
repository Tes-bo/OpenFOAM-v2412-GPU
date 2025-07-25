/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "faMesh.H"
#include "faMeshesRegistry.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshRegistry::faMeshRegistry
(
    const word& areaRegion,
    const polyMesh& mesh
)
:
    objectRegistry
    (
        IOobject
        (
            areaRegion,
            faMeshesRegistry::New(mesh).thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::fileName& Foam::faMeshRegistry::dbDir() const
{
    // In the usual case, the objectRegistry::dbDir() will be something
    // like finite-area/{region0,film} etc with the "finite-area/"
    // prefix coming from the enclosing registry of registries.
    //
    // So always check the name() portion, not the dbDir() itself
    // - either,    objectRegistry::dbDir().name()
    // - or (same), objectRegistry::name()

    if (objectRegistry::name() == polyMesh::defaultRegion)
    {
        return objectRegistry::parent().dbDir();
    }

    return objectRegistry::dbDir();
}


// ************************************************************************* //
