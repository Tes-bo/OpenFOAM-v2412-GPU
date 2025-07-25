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

#include "resolutionIndex.H"
#include "resolutionIndexModel.H"
#include "turbulenceModel.H"
#include "RASModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(resolutionIndex, 0);
    addToRunTimeSelectionTable(functionObject, resolutionIndex, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::resolutionIndex::resolutionIndex
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    resolutionIndexModelPtr_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::resolutionIndex::~resolutionIndex()
{}  // resolutionIndexModel was forward declared


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::resolutionIndex::read(const dictionary& dict)
{
    if (mesh_.foundObject<RASModelBase>(turbulenceModel::propertiesName))
    {
        FatalIOErrorInFunction(dict)
            << type() << " " << name()
            << " is not available for RANS-based turbulence models."
            << exit(FatalIOError);

        return false;
    }

    if (!fvMeshFunctionObject::read(dict))
    {
        return false;
    }

    Info<< type() << " " << name() << ":" << endl;

    resolutionIndexModelPtr_.reset
    (
        resolutionIndexModel::New(name(), mesh_, dict)
    );

    return true;
}


bool Foam::functionObjects::resolutionIndex::execute()
{
    if (!resolutionIndexModelPtr_->execute())
    {
        return false;
    }

    return true;
}


bool Foam::functionObjects::resolutionIndex::write()
{
    Info<< type() << " " << name() << " write:" << endl;

    if (!resolutionIndexModelPtr_->write())
    {
        return false;
    }
    Info<< endl;

    return true;
}


// ************************************************************************* //
