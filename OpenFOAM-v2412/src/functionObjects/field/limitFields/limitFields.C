/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd
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

#include "limitFields.H"
#include "fieldTypes.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(limitFields, 0);
    addToRunTimeSelectionTable(functionObject, limitFields, dictionary);
}
}

const Foam::Enum
<
    Foam::functionObjects::limitFields::limitType
>
Foam::functionObjects::limitFields::limitTypeNames_
({
    { limitType::CLAMP_MIN, "min" },
    { limitType::CLAMP_MAX, "max" },
    { limitType::CLAMP_RANGE, "range" },
    { limitType::CLAMP_RANGE, "both" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::functionObjects::limitFields::limitScalarField
(
    const word& fieldName
)
{
    auto* fieldPtr = obr_.getObjectPtr<volScalarField>(fieldName);
    if (!fieldPtr)
    {
        return false;
    }

    auto& field = *fieldPtr;

    Log << "    Limiting field " << fieldName << ":";

    if (limitType::CLAMP_NONE != withBounds_)
    {
        MinMax<scalar> currentRange = gMinMax(field);

        if (withBounds_ & limitType::CLAMP_MIN)
        {
            Log << " min(" << currentRange.min() << ')';
        }
        if (withBounds_ & limitType::CLAMP_MAX)
        {
            Log << " max(" << currentRange.max() << ')';
        }
    }

    if (limitType::CLAMP_MIN == withBounds_)
    {
        field.clamp_min(min_);
    }
    else if (limitType::CLAMP_MAX == withBounds_)
    {
        field.clamp_max(max_);
    }
    else if (limitType::CLAMP_RANGE == withBounds_)
    {
        field.clamp_range(min_, max_);
    }

    Log << endl;

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::limitFields::limitFields
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    fieldSet_(mesh_),
    withBounds_(limitType::CLAMP_NONE),
    min_(-VGREAT),
    max_(VGREAT)
{
    read(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::limitFields::read(const dictionary& dict)
{
    withBounds_ = limitType::CLAMP_NONE;

    if (fvMeshFunctionObject::read(dict))
    {
        Info<< type() << " " << name() << ":" << nl;

        withBounds_ = limitTypeNames_.get("limit", dict);

        if (withBounds_ & limitType::CLAMP_MIN)
        {
            min_ = dict.get<scalar>("min");
            Info<< "    Imposing lower limit " << min_ << nl;
        }

        if (withBounds_ & limitType::CLAMP_MAX)
        {
            max_ = dict.get<scalar>("max");
            Info<< "    Imposing upper limit " << max_ << nl;
        }

        fieldSet_.read(dict);

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::limitFields::execute()
{
    fieldSet_.updateSelection();

    Log << type() << ' ' << name() << ':' << nl;

    label count = 0, total = 0;

    for (const word& fieldName : fieldSet_.selectionNames())
    {
        ++total;
        if
        (
            limitScalarField(fieldName)
         || limitField<vector>(fieldName)
         || limitField<sphericalTensor>(fieldName)
         || limitField<symmTensor>(fieldName)
         || limitField<tensor>(fieldName)
        )
        {
            ++count;
        }
    }

    if (debug)
    {
        Log << " - limited " << count << '/' << total << " fields";
    }

    Log << endl;

    return true;
}


bool Foam::functionObjects::limitFields::write()
{
    for (const word& fieldName : fieldSet_.selectionNames())
    {
        lookupObject<regIOobject>(fieldName).write();
    }

    return true;
}


// ************************************************************************* //
