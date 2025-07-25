/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::multiply::initialiseResult(const word& fieldName)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    auto* fieldPtr = mesh_.cfindObject<volFieldType>(fieldName);

    if (fieldPtr)
    {
        auto* resultFieldPtr = mesh_.getObjectPtr<regIOobject>(resultName_);

        if (resultFieldPtr)
        {
            resultFieldPtr->checkOut();
        }

        Log << "    Initialising "
            << resultName_ << " to " << fieldPtr->name() << endl;

        return store(resultName_, tmp<volFieldType>::New(*fieldPtr));
    }

    return false;
}


template<class Type>
bool Foam::functionObjects::multiply::multiplyResult
(
    const word& fieldName,
    bool& processed
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

    auto* resultFieldPtr = mesh_.getObjectPtr<volFieldType>(resultName_);

    if (resultFieldPtr)
    {
        multiplyFieldType<Type, scalar>(*resultFieldPtr, fieldName, processed);
        multiplyFieldType<Type, vector>(*resultFieldPtr, fieldName, processed);
        multiplyFieldType<Type, sphericalTensor>
        (
            *resultFieldPtr,
            fieldName,
            processed
        );
        multiplyFieldType<Type, symmTensor>
        (
            *resultFieldPtr,
            fieldName,
            processed
        );
        multiplyFieldType<Type, tensor>(*resultFieldPtr, fieldName, processed);
    }

    return processed;
}


template<class Type1, class Type2>
typename std::enable_if
<
    Foam::functionObjects::multiply::is_valid_op<Type1, Type2>::value, bool
>::type
Foam::functionObjects::multiply::multiplyFieldType
(
    GeometricField<Type1, fvPatchField, volMesh>& result,
    const word& fieldName,
    bool& processed
)
{
    if (processed) return processed;

    typedef GeometricField<Type2, fvPatchField, volMesh> volFieldType;

    auto* fieldPtr = mesh_.cfindObject<volFieldType>(fieldName);

    if (fieldPtr)
    {
        Log << "    Performing " << result.name() << " * " << fieldPtr->name()
            << endl;

        auto newResult(result*(*fieldPtr));
        result.checkOut();

        store(resultName_, newResult);

        processed = true;
    }

    return processed;
}


template<class Type1, class Type2>
typename std::enable_if
<
    !Foam::functionObjects::multiply::is_valid_op<Type1, Type2>::value, bool
>::type
Foam::functionObjects::multiply::multiplyFieldType
(
    GeometricField<Type1, fvPatchField, volMesh>& result,
    const word& fieldName,
    bool& processed
)
{
    if (processed) return processed;

    typedef GeometricField<Type2, fvPatchField, volMesh> volFieldType;

    auto* fieldPtr = mesh_.cfindObject<volFieldType>(fieldName);

    if (fieldPtr)
    {
        Info<< "    Unsupported operation for "
            << result.name() << '(' << pTraits<Type1>::typeName << ')'
            << " * "
            << fieldPtr->name() << '(' << pTraits<Type2>::typeName << ')'
            << endl;
    }

    return processed;
}


// ************************************************************************* //
