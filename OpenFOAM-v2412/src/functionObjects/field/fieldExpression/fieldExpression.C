/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "fieldExpression.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(fieldExpression, 0);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::fieldExpression::setResultName
(
    const word& typeName,
    const word& defaultArg
)
{
    if (fieldName_.empty())
    {
        fieldName_ = defaultArg;
    }

    if (resultName_.empty())
    {
        if (fieldName_ != defaultArg)
        {
            resultName_ = typeName + '(' + fieldName_ + ')';
        }
        else
        {
            resultName_ = typeName;
        }

        resultName_ = scopedName(resultName_);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::fieldExpression::fieldExpression
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    const word& fieldName,
    const word& resultName
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldName_(fieldName),
    resultName_(resultName)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::fieldExpression::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        if (fieldName_.empty() || dict.found("field"))
        {
            dict.readEntry("field", fieldName_);
        }

        dict.readIfPresent("result", resultName_);

        if (dict.found("cellZones"))
        {
            zoneSubSetPtr_.reset(new Detail::zoneSubSet(mesh_, dict));
        }
        else
        {
            zoneSubSetPtr_.reset(nullptr);
        }

        return true;
    }

    return false;
}


bool Foam::functionObjects::fieldExpression::execute()
{
    if (!calc())
    {
        Warning
            << "    functionObjects::" << type() << " " << name()
            << " failed to execute." << endl;

        // Clear the result field from the objectRegistry if present
        clear();

        return false;
    }

    return true;
}


bool Foam::functionObjects::fieldExpression::write()
{
    return writeObject(resultName_);
}


bool Foam::functionObjects::fieldExpression::clear()
{
    return clearObject(resultName_);
}


// ************************************************************************* //
