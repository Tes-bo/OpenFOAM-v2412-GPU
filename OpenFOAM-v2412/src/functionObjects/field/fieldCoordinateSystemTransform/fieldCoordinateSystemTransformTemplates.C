/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2018 OpenCFD Ltd.
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

#include "fieldCoordinateSystemTransform.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transformGeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class FieldType>
void Foam::functionObjects::fieldCoordinateSystemTransform::transformField
(
    const FieldType& field
)
{
    word transFieldName(transformFieldName(field.name()));

    store
    (
        transFieldName,
        Foam::invTransform(dimensionedTensor(csysPtr_->R()), field)
    );
}


template<class FieldType, class RotationFieldType>
void Foam::functionObjects::fieldCoordinateSystemTransform::transformField
(
    const RotationFieldType& rot,
    const FieldType& field
)
{
    word transFieldName(transformFieldName(field.name()));

    store
    (
        transFieldName,
        Foam::invTransform(rot, field)
    );
}


template<class Type>
void Foam::functionObjects::fieldCoordinateSystemTransform::transform
(
    const word& fieldName
)
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;

    // Scalar quantities (bool, label, scalar) and sphericalTensor quantities
    // are transform invariant. Use (pTraits<Type>::nComponents == 1) to avoid
    // avoid generating a tensor field for a non-uniform transformation.

    if (foundObject<VolFieldType>(fieldName))
    {
        DebugInfo
            << type() << ": Field " << fieldName << " already in database"
            << endl;

        if (csysPtr_->uniform() || pTraits<Type>::nComponents == 1)
        {
            transformField<VolFieldType>
            (
                lookupObject<VolFieldType>(fieldName)
            );
        }
        else
        {
            transformField<VolFieldType>
            (
                vrotTensor(),
                lookupObject<VolFieldType>(fieldName)
            );
        }
    }
    else if (foundObject<SurfaceFieldType>(fieldName))
    {
        DebugInfo
            << type() << ": Field " << fieldName << " already in database"
            << endl;

        if (csysPtr_->uniform() || pTraits<Type>::nComponents == 1)
        {
            transformField<SurfaceFieldType>
            (
                lookupObject<SurfaceFieldType>(fieldName)
            );
        }
        else
        {
            transformField<SurfaceFieldType>
            (
                srotTensor(),
                lookupObject<SurfaceFieldType>(fieldName)
            );
        }
    }
    else
    {
        IOobject fieldHeader
        (
            fieldName,
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        );

        if (fieldHeader.typeHeaderOk<VolFieldType>(true, true, false))
        {
            DebugInfo
                << type() << ": Field " << fieldName << " read from file"
                << endl;

            if (csysPtr_->uniform() || pTraits<Type>::nComponents == 1)
            {
                transformField<VolFieldType>
                (
                    lookupObject<VolFieldType>(fieldName)
                );
            }
            else
            {
                transformField<VolFieldType>
                (
                    vrotTensor(),
                    lookupObject<VolFieldType>(fieldName)
                );
            }
        }
        else if (fieldHeader.typeHeaderOk<SurfaceFieldType>(true, true, false))
        {
            DebugInfo
                << type() << ": Field " << fieldName << " read from file"
                << endl;

            if (csysPtr_->uniform() || pTraits<Type>::nComponents == 1)
            {
                transformField<SurfaceFieldType>
                (
                    lookupObject<SurfaceFieldType>(fieldName)
                );
            }
            else
            {
                transformField<SurfaceFieldType>
                (
                    srotTensor(),
                    lookupObject<SurfaceFieldType>(fieldName)
                );
            }
        }
    }
}


// ************************************************************************* //
