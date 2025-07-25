/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

#include "cloud.H"
#include "Time.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(cloud, 0);
}

const Foam::word Foam::cloud::prefix("lagrangian");

Foam::word Foam::cloud::defaultName("defaultCloud");

const Foam::Enum<Foam::cloud::geometryType>
Foam::cloud::geometryTypeNames
({
    { geometryType::COORDINATES, "coordinates" },
    { geometryType::POSITIONS, "positions" }
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cloud::cloud(const objectRegistry& obr)
:
    cloud(obr, cloud::defaultName)
{}


Foam::cloud::cloud(const objectRegistry& obr, const word& cloudName)
:
    objectRegistry
    (
        IOobject
        (
            (cloudName.empty() ? cloud::defaultName : cloudName),
            obr.time().timeName(),
            cloud::prefix,
            obr,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            IOobject::REGISTER
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::cloud::nParcels() const
{
    NotImplemented;
    return 0;
}


void Foam::cloud::autoMap(const mapPolyMesh&)
{
    NotImplemented;
}


void Foam::cloud::readObjects(const objectRegistry& obr)
{
    NotImplemented;
}


void Foam::cloud::writeObjects(objectRegistry& obr) const
{
    NotImplemented;
}


// ************************************************************************* //
