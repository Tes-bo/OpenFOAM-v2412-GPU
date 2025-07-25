/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::limitFields::limitField(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;

    auto* fieldPtr = getObjectPtr<VolFieldType>(fieldName);
    if (!fieldPtr)
    {
        return false;
    }

    auto& field = *fieldPtr;

    Log << "    Limiting field " << fieldName << ":";

    const dimensionedScalar eps("eps", field.dimensions(), ROOTVSMALL);

    if (withBounds_ & limitType::CLAMP_MIN)
    {
        auto tmField = volScalarField::New
        (
            IOobject::scopedName(typeName, "mag" + field.name()),
            IOobject::NO_REGISTER,
            mag(field)
        );
        auto& mField = tmField.ref();

        Log << " min(|" << gMin(mField) << "|)";
        //field.normalise();
        field /= mag(field) + eps;
        mField.clamp_min(min_);
        field *= tmField;
    }

    if (withBounds_ & limitType::CLAMP_MAX)
    {
        auto tmField = volScalarField::New
        (
            IOobject::scopedName(typeName, "mag" + field.name()),
            IOobject::NO_REGISTER,
            mag(field)
        );
        auto& mField = tmField.ref();

        Log << " max(|" << gMax(mField) << "|)";
        //field.normalise();
        field /= mag(field) + eps;
        mField.clamp_max(max_);
        field *= tmField;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
