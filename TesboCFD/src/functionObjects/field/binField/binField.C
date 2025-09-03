/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "binField.H"
#include "binModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(binField, 0);
    addToRunTimeSelectionTable(functionObject, binField, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::binField::binField
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    binModelPtr_(nullptr)
{
    if (readFields)
    {
        read(dict);
    }
}


Foam::functionObjects::binField::binField
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    binModelPtr_(nullptr)
{
    if (readFields)
    {
        read(dict);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::binField::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    Info<< type() << " " << name() << ":" << endl;

    binModelPtr_.reset(binModel::New(dict, mesh_, name()));

    return true;
}


bool Foam::functionObjects::binField::execute()
{
    Log << type() << " " << name() << ":" << nl
        << "    Calculating bins" << nl << endl;

    binModelPtr_->apply();

    return true;
}


bool Foam::functionObjects::binField::write()
{
    return true;
}


void Foam::functionObjects::binField::updateMesh(const mapPolyMesh& mpm)
{
    binModelPtr_->updateMesh(mpm);
}


void Foam::functionObjects::binField::movePoints(const polyMesh& mesh)
{
    binModelPtr_->movePoints(mesh);
}


// ************************************************************************* //
