/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2019 OpenCFD Ltd.
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

#include "fieldAverageItem.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "polySurfaceFields.H"
#include "OFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::fieldAverage::addMeanFieldType
(
    fieldAverageItem& item
)
{
    const Type* fieldPtr = findObject<Type>(item.fieldName());

    if (!fieldPtr)
    {
        return false;
    }

    // Field has been found, so set active flag to true
    item.active() = true;

    const word& meanFieldName = item.meanFieldName();

    Log << "    Reading/initialising field " << meanFieldName << endl;

    if (foundObject<Type>(meanFieldName))
    {}
    else if (obr().found(meanFieldName))
    {
        Log << "    Cannot allocate average field " << meanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        item.mean() = false;
    }
    else
    {
        const Type& baseField = *fieldPtr;

        // Store on registry
        obr().store
        (
            new Type
            (
                IOobject
                (
                    meanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    (
                        restartOnOutput_
                      ? IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT
                    ),
                    IOobject::NO_WRITE
                ),
                1*baseField
            )
        );

        return true;
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::fieldAverage::addMeanField
(
    fieldAverageItem& item
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal VolFieldInternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    bool added = false;

    if (item.mean())
    {
        added =
        (
            addMeanFieldType<VolFieldType>(item)
         || addMeanFieldType<VolFieldInternalType>(item)
         || addMeanFieldType<SurfaceFieldType>(item)
         || addMeanFieldType<SurfFieldType>(item)
        );
    }

    return added;
}


template<class Type>
bool Foam::functionObjects::fieldAverage::restoreWindowFieldsType
(
    const fieldAverageItem& item
)
{
    if (restartOnOutput_)
    {
        return false;
    }

    const word& fieldName = item.fieldName();

    const Type* fieldPtr = findObject<Type>(fieldName);

    if (!fieldPtr)
    {
        return false;
    }

    const FIFOStack<word>& fieldNames = item.windowFieldNames();

    forAllConstIters(fieldNames, fieldIter)
    {
        const word& name = fieldIter();

        IOobject io
        (
            name,
            obr().time().timeName(obr().time().startTime().value()),
            obr(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::REGISTER
        );

        if (io.typeHeaderOk<Type>(true))
        {
            DebugInfo << "Read and store: " << name << endl;
            obr().store(new Type(io, fieldPtr->mesh()));
        }
        else
        {
            WarningInFunction
                << "Unable to read window " << Type::typeName << " " << name
                << ".  Averaging restart behaviour may be compromised"
                << endl;
        }
    }

    return true;
}


template<class Type>
void Foam::functionObjects::fieldAverage::restoreWindowFields
(
    const fieldAverageItem& item
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal VolFieldInternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    if (item.window() > 0)
    {
        (void)
        (
            restoreWindowFieldsType<VolFieldType>(item)
         || restoreWindowFieldsType<VolFieldInternalType>(item)
         || restoreWindowFieldsType<SurfaceFieldType>(item)
         || restoreWindowFieldsType<SurfFieldType>(item)
        );
    }
}


template<class Type1, class Type2>
bool Foam::functionObjects::fieldAverage::addPrime2MeanFieldType
(
    fieldAverageItem& item
)
{
    const auto* baseFieldPtr = findObject<Type1>(item.fieldName());

    if (!baseFieldPtr)
    {
        return false;
    }

    const word& meanFieldName = item.meanFieldName();
    const word& prime2MeanFieldName = item.prime2MeanFieldName();

    Log << "    Reading/initialising field " << prime2MeanFieldName << nl;

    if (foundObject<Type2>(prime2MeanFieldName))
    {}
    else if (obr().found(prime2MeanFieldName))
    {
        Log << "    Cannot allocate average field " << prime2MeanFieldName
            << " since an object with that name already exists."
            << " Disabling averaging for field." << endl;

        item.prime2Mean() = false;
    }
    else
    {
        const auto& baseField = *baseFieldPtr;
        const Type1& meanField = lookupObject<Type1>(meanFieldName);

        // Store on registry
        obr().store
        (
            new Type2
            (
                IOobject
                (
                    prime2MeanFieldName,
                    obr().time().timeName(obr().time().startTime().value()),
                    obr(),
                    restartOnOutput_?
                        IOobject::NO_READ
                      : IOobject::READ_IF_PRESENT,
                    IOobject::NO_WRITE
                ),
                sqr(baseField) - sqr(meanField)
            )
        );

        return true;
    }

    return false;
}


template<class Type1, class Type2>
bool Foam::functionObjects::fieldAverage::addPrime2MeanField
(
    fieldAverageItem& item
)
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal VolFieldInternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal VolFieldInternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    bool added = false;

    if (item.prime2Mean())
    {
        if (!item.mean())
        {
            FatalErrorInFunction
                << "To calculate the prime-squared average, the "
                << "mean average must also be selected for field "
                << item.fieldName() << nl << exit(FatalError);
        }

        added =
            addPrime2MeanFieldType<VolFieldType1, VolFieldType2>(item)
         || addPrime2MeanFieldType<VolFieldInternalType1, VolFieldInternalType2>
            (
                item
            )
         || addPrime2MeanFieldType<SurfaceFieldType1, SurfaceFieldType2>(item)
         || addPrime2MeanFieldType<SurfFieldType1, SurfFieldType2>(item);
    }

    return added;
}


template<class Type>
bool Foam::functionObjects::fieldAverage::storeWindowFieldType
(
    fieldAverageItem& item
)
{
    const auto* fPtr = findObject<Type>(item.fieldName());

    if (!fPtr)
    {
        return false;
    }

    const Type& baseField = *fPtr;

    const word windowFieldName = item.windowFieldName(this->name());

    // Store on registry
    obr().store
    (
        new Type
        (
            IOobject
            (
                windowFieldName,
                obr().time().timeName(obr().time().startTime().value()),
                obr(),
                restartOnOutput_ ?
                    IOobject::NO_READ
                  : IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            1*baseField
        )
    );

    DebugInfo << "Create and store: " << windowFieldName << endl;

    item.addToWindow(windowFieldName, obr().time().deltaTValue());

    return true;
}


template<class Type>
void Foam::functionObjects::fieldAverage::storeWindowFields()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal VolFieldInternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    for (fieldAverageItem& item : faItems_)
    {
        if (item.storeWindowFields())
        {
            (void)
            (
                storeWindowFieldType<VolFieldType>(item)
             || storeWindowFieldType<VolFieldInternalType>(item)
             || storeWindowFieldType<SurfaceFieldType>(item)
             || storeWindowFieldType<SurfFieldType>(item)
            );
        }
    }
}


template<class Type>
void Foam::functionObjects::fieldAverage::calculateMeanFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal VolFieldInternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    const auto& obr = this->obr();

    for (const fieldAverageItem& item : faItems_)
    {
        (void)
        (
            item.calculateMeanField<VolFieldType>(obr)
         || item.calculateMeanField<VolFieldInternalType>(obr)
         || item.calculateMeanField<SurfaceFieldType>(obr)
         || item.calculateMeanField<SurfFieldType>(obr)
        );
    }
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::calculatePrime2MeanFields() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal VolFieldInternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal VolFieldInternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    const auto& obr = this->obr();

    for (const fieldAverageItem& item : faItems_)
    {
        (void)
        (
            item.calculatePrime2MeanField<VolFieldType1, VolFieldType2>(obr)
         || item.calculatePrime2MeanField
            <VolFieldInternalType1, VolFieldInternalType2>(obr)
         || item.calculatePrime2MeanField
            <SurfaceFieldType1, SurfaceFieldType2>(obr)
         || item.calculatePrime2MeanField<SurfFieldType1, SurfFieldType2>(obr)
        );
    }
}


template<class Type1, class Type2>
bool Foam::functionObjects::fieldAverage::addMeanSqrToPrime2MeanType
(
    const fieldAverageItem& item
) const
{
    if (!foundObject<Type1>(item.fieldName()))
    {
        return false;
    }

    const Type1& meanField = lookupObject<Type1>(item.meanFieldName());

    Type2& prime2MeanField = lookupObjectRef<Type2>(item.prime2MeanFieldName());

    prime2MeanField += sqr(meanField);

    return true;
}


template<class Type1, class Type2>
void Foam::functionObjects::fieldAverage::addMeanSqrToPrime2Mean() const
{
    typedef GeometricField<Type1, fvPatchField, volMesh> VolFieldType1;
    typedef typename VolFieldType1::Internal VolFieldInternalType1;
    typedef GeometricField<Type1, fvsPatchField, surfaceMesh> SurfaceFieldType1;
    typedef DimensionedField<Type1, polySurfaceGeoMesh> SurfFieldType1;

    typedef GeometricField<Type2, fvPatchField, volMesh> VolFieldType2;
    typedef typename VolFieldType2::Internal VolFieldInternalType2;
    typedef GeometricField<Type2, fvsPatchField, surfaceMesh> SurfaceFieldType2;
    typedef DimensionedField<Type2, polySurfaceGeoMesh> SurfFieldType2;

    for (const fieldAverageItem& item : faItems_)
    {
        if (item.prime2Mean())
        {
            (void)
            (
                addMeanSqrToPrime2MeanType<VolFieldType1, VolFieldType2>(item)
             || addMeanSqrToPrime2MeanType
                <VolFieldInternalType1, VolFieldInternalType2>(item)
             || addMeanSqrToPrime2MeanType
                <SurfaceFieldType1, SurfaceFieldType2>(item)
             || addMeanSqrToPrime2MeanType
                <SurfFieldType1, SurfFieldType2>(item)
            );
        }
    }
}


template<class Type>
bool Foam::functionObjects::fieldAverage::writeFieldType
(
    const word& fieldName
) const
{
    const auto* fPtr = findObject<Type>(fieldName);

    if (fPtr)
    {
        DebugInfo<< "writing " << Type::typeName << ": " << fieldName << endl;
        return fPtr->write();
    }

    return false;
}


template<class Type>
void Foam::functionObjects::fieldAverage::writeFields() const
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef typename VolFieldType::Internal VolFieldInternalType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    for (const fieldAverageItem& item : faItems_)
    {
        if (item.mean())
        {
            const word& fieldName = item.meanFieldName();

            (void)
            (
                writeFieldType<VolFieldType>(fieldName)
             || writeFieldType<VolFieldInternalType>(fieldName)
             || writeFieldType<SurfaceFieldType>(fieldName)
             || writeFieldType<SurfFieldType>(fieldName)
            );
        }

        if (item.prime2Mean())
        {
            const word& fieldName = item.prime2MeanFieldName();

            (void)
            (
                writeFieldType<VolFieldType>(fieldName)
             || writeFieldType<VolFieldInternalType>(fieldName)
             || writeFieldType<SurfaceFieldType>(fieldName)
             || writeFieldType<SurfFieldType>(fieldName)
            );
        }

        if (item.writeWindowFields())
        {
            FIFOStack<word> fieldNames = item.windowFieldNames();
            forAllConstIters(fieldNames, fieldNameIter)
            {
                const word& fieldName = fieldNameIter();

                (void)
                (
                    writeFieldType<VolFieldType>(fieldName)
                 || writeFieldType<VolFieldInternalType>(fieldName)
                 || writeFieldType<SurfaceFieldType>(fieldName)
                 || writeFieldType<SurfFieldType>(fieldName)
                );
            }
        }
    }
}


// ************************************************************************* //
