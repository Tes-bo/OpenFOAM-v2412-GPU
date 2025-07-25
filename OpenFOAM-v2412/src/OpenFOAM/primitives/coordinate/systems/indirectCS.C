/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "indirectCS.H"
#include "coordinateSystems.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace coordSystem
{
    defineTypeNameAndDebug(indirect, 0);
    addToRunTimeSelectionTable(coordinateSystem, indirect, registry);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::coordSystem::indirect::indirect(const indirect& csys)
:
    coordinateSystem(),
    backend_(csys.backend_)
{}


Foam::coordSystem::indirect::indirect(indirect&& csys)
:
    coordinateSystem(),
    backend_(std::move(csys.backend_))
{}


// Use lookup() instead of cfind() to trigger FatalError on any problems
Foam::coordSystem::indirect::indirect
(
    const objectRegistry& obr,
    const word& name
)
:
    coordinateSystem(),
    backend_(&(coordinateSystems::New(obr).lookup(name)))
{}


Foam::coordSystem::indirect::indirect
(
    const objectRegistry& obr,
    const dictionary& dict,
    IOobjectOption::readOption /* (unused) */
)
:
    indirect(obr, dict.get<word>("name"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::coordSystem::indirect::write(Ostream& os) const
{
    writeEntry(coordinateSystem::typeName, os);
}


void Foam::coordSystem::indirect::writeEntry(Ostream& os) const
{
    writeEntry(coordinateSystem::typeName, os);
}


void Foam::coordSystem::indirect::writeEntry
(
    const word& keyword,
    Ostream& os
) const
{
    if (!good())
    {
        return;
    }

    const bool subDict = !keyword.empty();

    if (subDict)
    {
        os.beginBlock(keyword);

        os.writeEntry("type", type());
        os.writeEntry("name", name());

        os.endBlock();
    }
}


// ************************************************************************* //
