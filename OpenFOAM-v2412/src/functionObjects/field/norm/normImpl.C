/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::functionObjects::norm::calcNorm()
{
    typedef GeometricField<Type, fvPatchField, volMesh> VolFieldType;
    typedef GeometricField<Type, fvsPatchField, surfaceMesh> SurfaceFieldType;
    typedef DimensionedField<Type, polySurfaceGeoMesh> SurfFieldType;

    if (mesh_.foundObject<VolFieldType>(fieldName_))
    {
        return store
        (
            resultName_,
            calcNormType<VolFieldType>()
        );
    }
    else if (mesh_.foundObject<SurfaceFieldType>(fieldName_))
    {
        return store
        (
            resultName_,
            calcNormType<SurfaceFieldType>()
        );
    }
    else if (mesh_.foundObject<SurfFieldType>(fieldName_))
    {
        return store
        (
            resultName_,
            calcNormType<SurfFieldType>()
        );
    }

    return false;
}


template<class GeoFieldType>
Foam::tmp<GeoFieldType> Foam::functionObjects::norm::calcNormType()
{
    const GeoFieldType& field = mesh_.lookupObject<GeoFieldType>(fieldName_);

    const dimensionedScalar perturb(field.dimensions(), SMALL);

    switch (norm_)
    {
        case normType::L1:
        {
            return field/stabilise(sumMag(field), perturb);
        }

        case normType::L2:
        {
            return field/stabilise(mag(field), perturb);
        }

        case normType::LP:
        {
            return
                field
               /stabilise
                (
                   pow(pow(mag(field), p_), scalar(1)/p_),
                   perturb
                );
        }

        case normType::MAX:
        {
            return field/stabilise(max(mag(field)), perturb);
        }

        case normType::COMPOSITE:
        {
            const scalar t = mesh_.time().timeOutputValue();

            const dimensionedScalar divisor
            (
                field.dimensions(),
                divisorPtr_->value(t)
            );

            return field/stabilise(divisor, perturb);
        }

        case normType::FIELD:
        {
            return field/stabilise(fieldNorm(field), perturb);
        }

        default:
            break;
    }

    return nullptr;
}


template<class Type>
Foam::tmp<Foam::volScalarField> Foam::functionObjects::norm::fieldNorm
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    return mesh_.lookupObject<volScalarField>(divisorFieldName_);
}


template<class Type>
Foam::tmp<Foam::surfaceScalarField> Foam::functionObjects::norm::fieldNorm
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>&
)
{
    return mesh_.lookupObject<surfaceScalarField>(divisorFieldName_);
}


template<class Type>
Foam::tmp<Foam::polySurfaceScalarField> Foam::functionObjects::norm::fieldNorm
(
    const DimensionedField<Type, polySurfaceGeoMesh>&
)
{
    return mesh_.lookupObject<polySurfaceScalarField>(divisorFieldName_);
}


// ************************************************************************* //
