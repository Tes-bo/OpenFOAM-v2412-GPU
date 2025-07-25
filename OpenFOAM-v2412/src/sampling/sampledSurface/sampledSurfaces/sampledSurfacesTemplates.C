/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2024 OpenCFD Ltd.
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

#include "sampledSurfaces.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "polySurfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::sampledSurfaces::writeSurface
(
    surfaceWriter& writer,
    const Field<Type>& values,
    const word& fieldName
)
{
    fileName outputName = writer.write(fieldName, values);

    // Case-local file name with "<case>" to make relocatable

    dictionary propsDict;
    propsDict.add
    (
        "file",
        time_.relativePath(outputName, true)
    );
    setProperty(fieldName, propsDict);
}


template<class Type, class GeoMeshType>
void Foam::sampledSurfaces::storeRegistryField
(
    const sampledSurface& s,
    const word& fieldName,
    const dimensionSet& dims,
    Field<Type>&& values
)
{
    s.storeRegistryField<Type, GeoMeshType>
    (
        storedObjects(),
        fieldName,
        dims,
        std::move(values),
        IOobject::groupName(name(), s.name())  // surfaceName
    );
}


template<class Type>
void Foam::sampledSurfaces::performAction
(
    const VolumeField<Type>& fld,
    unsigned request
)
{
    // The sampler for this field
    autoPtr<interpolation<Type>> samplePtr;

    // The interpolator for this field
    autoPtr<interpolation<Type>> interpPtr;

    const word& fieldName = fld.name();

    const dimensionSet& dims = fld.dimensions();

    forAll(*this, surfi)
    {
        // Skip empty surfaces (eg, failed cut-plane)
        if (!hasFaces_[surfi]) continue;

        const sampledSurface& s = operator[](surfi);

        Field<Type> values;

        if (s.isPointData())
        {
            if (!interpPtr)
            {
                interpPtr = interpolation<Type>::New
                (
                    sampleNodeScheme_,
                    fld
                );
            }

            values = s.interpolate(*interpPtr);
        }
        else
        {
            if (!samplePtr)
            {
                samplePtr = interpolation<Type>::New
                (
                    sampleFaceScheme_,
                    fld
                );
            }

            values = s.sample(*samplePtr);
        }

        if ((request & actions_[surfi]) & ACTION_WRITE)
        {
            writeSurface<Type>(writers_[surfi], values, fieldName);
        }

        if ((request & actions_[surfi]) & ACTION_STORE)
        {
            if (s.isPointData())
            {
                storeRegistryField<Type, polySurfacePointGeoMesh>
                (
                    s, fieldName, dims, std::move(values)
                );
            }
            else
            {
                storeRegistryField<Type, polySurfaceGeoMesh>
                (
                    s, fieldName, dims, std::move(values)
                );
            }
        }
    }
}


template<class Type>
void Foam::sampledSurfaces::performAction
(
    const SurfaceField<Type>& fld,
    unsigned request
)
{
    const word& fieldName = fld.name();

    const dimensionSet& dims = fld.dimensions();

    forAll(*this, surfi)
    {
        // Skip empty surfaces (eg, failed cut-plane)
        if (!hasFaces_[surfi]) continue;

        const sampledSurface& s = (*this)[surfi];

        Field<Type> values(s.sample(fld));

        if ((request & actions_[surfi]) & ACTION_WRITE)
        {
            writeSurface<Type>(writers_[surfi], values, fieldName);
        }

        if ((request & actions_[surfi]) & ACTION_STORE)
        {
            storeRegistryField<Type, polySurfaceGeoMesh>
            (
                s, fieldName, dims, std::move(values)
            );
        }
    }
}


template<class GeoField>
void Foam::sampledSurfaces::performAction
(
    const IOobjectList& objects,
    unsigned request
)
{
    wordList fieldNames;
    if (loadFromFiles_)
    {
        fieldNames = objects.sortedNames<GeoField>(fieldSelection_);
    }
    else
    {
        fieldNames = mesh_.thisDb().sortedNames<GeoField>(fieldSelection_);
    }

    for (const word& fieldName : fieldNames)
    {
        if (verbose_)
        {
            Info<< "sampleWrite: " << fieldName << endl;
        }

        refPtr<GeoField> tfield;

        if (loadFromFiles_)
        {
            tfield.emplace
            (
                IOobject
                (
                    fieldName,
                    time_.timeName(),
                    mesh_.thisDb(),
                    IOobjectOption::MUST_READ,
                    IOobjectOption::NO_WRITE,
                    IOobjectOption::NO_REGISTER
                ),
                mesh_
            );
        }
        else
        {
            tfield.cref(mesh_.thisDb().cfindObject<GeoField>(fieldName));
        }

        if (tfield)
        {
            performAction(tfield(), request);
        }
    }
}


template<class Container, class Predicate>
bool Foam::sampledSurfaces::testAny
(
    const Container& items,
    const Predicate& pred
)
{
    for (const auto& item : items)
    {
        if (pred(item))
        {
            return true;
        }
    }

    return false;
}


// ************************************************************************* //
