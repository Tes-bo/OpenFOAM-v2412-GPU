/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "fieldValue.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldValue, 0);
    defineRunTimeSelectionTable(fieldValue, runTime);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& valueType
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, valueType, dict),
    writeFields_(false),
    regionName_(),
    scaleFactor_(1.0),
    dict_(dict),
    fields_()
{
    read(dict);
}


Foam::functionObjects::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const word& valueType
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(obr_, name, valueType, dict),
    writeFields_(false),
    regionName_(),
    scaleFactor_(1.0),
    dict_(dict),
    fields_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldValue::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict) && writeFile::read(dict))
    {
        if (dict != dict_)
        {
            dict_ = dict;
        }

        dict.readEntry("writeFields", writeFields_);
        scaleFactor_ = dict.getOrDefault<scalar>("scaleFactor", 1.0);
        dict.readEntry("fields", fields_);

        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldValue::execute()
{
    return true;
}


bool Foam::functionObjects::fieldValue::write()
{
    Log << type() << ' ' << name() << " write:" << nl;

    return true;
}


// ************************************************************************* //
