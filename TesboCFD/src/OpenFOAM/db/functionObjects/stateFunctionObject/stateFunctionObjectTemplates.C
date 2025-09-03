/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Type Foam::functionObjects::stateFunctionObject::getProperty
(
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getProperty(entryName, result);
    return result;
}


template<class Type>
bool Foam::functionObjects::stateFunctionObject::getProperty
(
    const word& entryName,
    Type& value
) const
{
    return getObjectProperty(name(), entryName, value);
}


template<class Type>
void Foam::functionObjects::stateFunctionObject::setProperty
(
    const word& entryName,
    const Type& value
)
{
    setObjectProperty(name(), entryName, value);
}


template<class Type>
Type Foam::functionObjects::stateFunctionObject::getObjectProperty
(
    const word& objectName,
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    getObjectProperty(objectName, entryName, result);
    return result;
}


template<class Type>
bool Foam::functionObjects::stateFunctionObject::getObjectProperty
(
    const word& objectName,
    const word& entryName,
    Type& value
) const
{
    return stateDict().getObjectProperty(objectName, entryName, value);
}


template<class Type>
void Foam::functionObjects::stateFunctionObject::setObjectProperty
(
    const word& objectName,
    const word& entryName,
    const Type& value
)
{
    stateDict().setObjectProperty(objectName, entryName, value);
}


template<class Type>
void Foam::functionObjects::stateFunctionObject::setResult
(
    const word& entryName,
    const Type& value
)
{
    setObjectResult(name(), entryName, value);
}


template<class Type>
void Foam::functionObjects::stateFunctionObject::setObjectResult
(
    const word& objectName,
    const word& entryName,
    const Type& value
)
{
    stateDict().setObjectResult(objectName, entryName, value);
}


template<class Type>
Type Foam::functionObjects::stateFunctionObject::getResult
(
    const word& entryName,
    const Type& defaultValue
) const
{
    return getObjectResult(name(), entryName, defaultValue);
}


template<class Type>
Type Foam::functionObjects::stateFunctionObject::getObjectResult
(
    const word& objectName,
    const word& entryName,
    const Type& defaultValue
) const
{
    Type result = defaultValue;
    (void)getObjectResult(objectName, entryName, result);
    return result;
}


template<class Type>
bool Foam::functionObjects::stateFunctionObject::getObjectResult
(
    const word& objectName,
    const word& entryName,
    Type& value
) const
{
    return stateDict().getObjectResult(objectName, entryName, value);
}


// ************************************************************************* //
