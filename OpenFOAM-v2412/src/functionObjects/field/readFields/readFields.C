/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "readFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(readFields, 0);
    addToRunTimeSelectionTable(functionObject, readFields, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::readFields::readFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    readOnStart_(dict.getOrDefault("readOnStart", true)),
    fieldSet_()
{
    read(dict);

    if (readOnStart_)
    {
        execute();
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::readFields::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    dict.readEntry("fields", fieldSet_);

    return true;
}


bool Foam::functionObjects::readFields::execute()
{
    for (const word& fieldName : fieldSet_)
    {
        // Already loaded?
        auto* ptr = mesh_.getObjectPtr<regIOobject>(fieldName);

        if (ptr)
        {
            if (functionObject::postProcess)
            {
                DebugInfo
                    << "readFields : "
                    << ptr->name() << " (" << ptr->type()
                    << ") already in database - removing" << endl;

                ptr->checkOut();
            }
            else
            {
                DebugInfo
                    << "readFields : "
                    << ptr->name() << " (" << ptr->type()
                    << ") already in database" << endl;

                continue;
            }
        }

        // Load field as necessary
        IOobject io
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        );

        const bool ok =
        (
            io.typeHeaderOk<regIOobject>(false)
         &&
            (
                loadField<scalar>(io)
             || loadField<vector>(io)
             || loadField<sphericalTensor>(io)
             || loadField<symmTensor>(io)
             || loadField<tensor>(io)
            )
        );

        if (!ok)
        {
            DebugInfo
                << "readFields : failed to load " << fieldName << endl;
        }
    }

    return true;
}


bool Foam::functionObjects::readFields::write()
{
    return true;
}


// ************************************************************************* //
