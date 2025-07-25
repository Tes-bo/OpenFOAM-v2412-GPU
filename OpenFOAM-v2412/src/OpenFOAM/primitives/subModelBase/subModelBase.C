/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "subModelBase.H"

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

bool Foam::subModelBase::subModelBase::inLine() const
{
    return (!modelName_.empty());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::subModelBase::subModelBase(dictionary& properties)
:
    modelName_(),
    properties_(properties),
    dict_(),
    baseName_(),
    modelType_(),
    coeffDict_(),
    log(properties.getOrDefault<bool>("log", true))
{}


Foam::subModelBase::subModelBase
(
    dictionary& properties,
    const dictionary& dict,
    const word& baseName,
    const word& modelType,
    const word& dictExt
)
:
    modelName_(),
    properties_(properties),
    dict_(dict),
    baseName_(baseName),
    modelType_(modelType),
    coeffDict_(dict.subDict(modelType + dictExt)),
    log(coeffDict_.getOrDefault<bool>("log", true))
{}


Foam::subModelBase::subModelBase
(
    const word& modelName,
    dictionary& properties,
    const dictionary& dict,
    const word& baseName,
    const word& modelType
)
:
    modelName_(modelName),
    properties_(properties),
    dict_(dict),
    baseName_(baseName),
    modelType_(modelType),
    coeffDict_(dict),
    log(coeffDict_.getOrDefault<bool>("log", true))
{}


Foam::subModelBase::subModelBase(const subModelBase& smb)
:
    modelName_(smb.modelName_),
    properties_(smb.properties_),
    dict_(smb.dict_),
    baseName_(smb.baseName_),
    modelType_(smb.modelType_),
    coeffDict_(smb.coeffDict_),
    log(coeffDict_.getOrDefault<bool>("log", true))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::word& Foam::subModelBase::modelName() const
{
    return modelName_;
}


const Foam::dictionary& Foam::subModelBase::dict() const
{
    return dict_;
}


const Foam::word& Foam::subModelBase::baseName() const
{
    return baseName_;
}


const Foam::word& Foam::subModelBase::modelType() const
{
    return modelType_;
}


const Foam::dictionary& Foam::subModelBase::coeffDict() const
{
    return coeffDict_;
}


const Foam::dictionary& Foam::subModelBase::properties() const
{
    return properties_;
}


bool Foam::subModelBase::defaultCoeffs(const bool printMsg) const
{
    bool def = coeffDict_.getOrDefault("defaultCoeffs", false);
    if (printMsg && def)
    {
        // Note: not using Log<< for output
        Info<< incrIndent;
        Info<< indent << "Employing default coefficients" << endl;
        Info<< decrIndent;
    }

    return def;
}


bool Foam::subModelBase::active() const
{
    return true;
}


void Foam::subModelBase::cacheFields(const bool)
{}


bool Foam::subModelBase::writeTime() const
{
    return active();
}


Foam::fileName Foam::subModelBase::localPath() const
{
    if (!modelName_.empty())
    {
        return modelName_;
    }

    return baseName_;
}


bool Foam::subModelBase::getModelDict
(
    const word& entryName,
    dictionary& dict
) const
{
    if (properties_.found(baseName_))
    {
        const dictionary& baseDict = properties_.subDict(baseName_);

        if (inLine() && baseDict.found(modelName_))
        {
            const dictionary& modelDict = baseDict.subDict(modelName_);
            dict = modelDict.subOrEmptyDict(entryName);
            return true;
        }
        else if (baseDict.found(modelType_))
        {
            const dictionary& modelDict = baseDict.subDict(modelType_);
            dict = modelDict.subOrEmptyDict(entryName);
            return true;
        }
    }

    return false;
}


void Foam::subModelBase::write(Ostream& os) const
{
    os  << coeffDict_;
}


// ************************************************************************* //
