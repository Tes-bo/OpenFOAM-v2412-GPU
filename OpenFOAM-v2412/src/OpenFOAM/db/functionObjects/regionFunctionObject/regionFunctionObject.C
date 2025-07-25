/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "regionFunctionObject.H"
#include "Time.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(regionFunctionObject, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

const Foam::objectRegistry&
Foam::functionObjects::regionFunctionObject::obr() const
{
    if (!obrPtr_ && !subRegistryName_.empty())
    {
        // Recursive - so we also find things registered on Time
        obrPtr_ = obr_.cfindObject<objectRegistry>(subRegistryName_, true);

        // Also search functionObject output ("functionObjectObjects")
        if (!obrPtr_)
        {
            obrPtr_ =
                storedObjects().cfindObject<objectRegistry>(subRegistryName_);
        }

        if (!obrPtr_)
        {
            WarningInFunction
                << "Could not locate subRegion \""
                << subRegistryName_ << '"' << nl;
        }
    }

    return (obrPtr_ ? *obrPtr_ : obr_);
}


bool Foam::functionObjects::regionFunctionObject::writeObject
(
    const word& fieldName
)
{
    const regIOobject* ptr = this->cfindObject<regIOobject>(fieldName);

    if (ptr)
    {
        Log << "    functionObjects::" << type() << " " << name()
            << " writing field: " << ptr->name() << endl;

        ptr->write();

        return true;
    }

    return false;
}


bool Foam::functionObjects::regionFunctionObject::clearObject
(
    const word& fieldName
)
{
    // Same as getObjectPtr, since the object is already non-const
    regIOobject* ptr = this->findObject<regIOobject>(fieldName);

    if (ptr)
    {
        if (ptr->ownedByRegistry())
        {
            return ptr->checkOut();
        }
        else
        {
            return false;
        }
    }

    return true;
}


void Foam::functionObjects::regionFunctionObject::clearObjects
(
    const wordList& objNames
)
{
    for (const word& objName : objNames)
    {
        regIOobject* ptr = this->findObject<regIOobject>(objName);

        if (ptr && ptr->ownedByRegistry())
        {
            ptr->checkOut();
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    stateFunctionObject(name, runTime),
    subRegistryName_(dict.getOrDefault<word>("subRegion", word::null)),
    obr_
    (
        runTime.lookupObject<objectRegistry>
        (
            dict.getOrDefault("region", polyMesh::defaultRegion)
        )
    ),
    obrPtr_(nullptr)
{}


Foam::functionObjects::regionFunctionObject::regionFunctionObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    stateFunctionObject(name, obr.time()),
    subRegistryName_(dict.getOrDefault<word>("subRegion", word::null)),
    obr_(obr),
    obrPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::regionFunctionObject::read(const dictionary& dict)
{
    stateFunctionObject::read(dict);

    subRegistryName_ = dict.getOrDefault<word>("subRegion", word::null);

    obrPtr_ = nullptr;

    return true;
}


// ************************************************************************* //
