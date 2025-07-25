/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const IOobject* io,
    const typename GeoField::Mesh& mesh,
    const bool syncPar
)
{
    if (io)
    {
        return tmp<GeoField>::New(*io, mesh);
    }

    return nullptr;
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const IOobject* io,
    const fvMeshSubsetProxy& proxy,
    const bool syncPar,
    objectRegistry* cache
)
{
    tmp<GeoField> tfield;

    if (io)
    {
        const word& fieldName = io->name();

        if (cache)
        {
            // Get reference from cache if possible
            tfield.cref(cache->cfindObject<GeoField>(fieldName));

            if (tfield)
            {
                return tfield;
            }
        }

        tfield = proxy.interpolate
        (
            getField<GeoField>(io, proxy.baseMesh(), syncPar)
        );

        if (tfield && cache)
        {
            // Move field to the cache
            IOobject newIO(tfield(), *cache);
            newIO.readOpt(IOobjectOption::NO_READ);
            newIO.writeOpt(IOobjectOption::NO_WRITE);

            tfield.ref().checkOut();  // Paranoid
            cache->store(new GeoField(newIO, tfield));

            tfield.cref(cache->cfindObject<GeoField>(fieldName));
        }
    }

    return tfield;
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects,
    const word& fieldName,
    const bool syncPar
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>
    (
        objects.findObject(fieldName),
        mesh,
        syncPar
    );
}


template<class GeoField>
Foam::tmp<GeoField> Foam::getField
(
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects,
    const word& fieldName,
    const bool syncPar,
    objectRegistry* cache
)
{
    // Can do something with syncPar on failure ...

    return getField<GeoField>
    (
        objects.findObject(fieldName),
        proxy,
        syncPar,
        cache
    );
}


template<class GeoField>
Foam::PtrList<const GeoField> Foam::readFields
(
    const typename GeoField::Mesh& mesh,
    const IOobjectList& objects
)
{
    const bool syncPar = true;

    // Available fields of type GeoField, sorted order
    const wordList fieldNames(objects.sortedNames<GeoField>());

    // Construct the fields
    PtrList<const GeoField> fields(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        auto tfield =
            getField<GeoField>(mesh, objects, fieldName, syncPar);

        if (tfield)
        {
            fields.set(nFields++, tfield.ptr());
        }
    }

    fields.resize(nFields);
    return fields;
}


template<class GeoField>
Foam::PtrList<const GeoField> Foam::readFields
(
    const fvMeshSubsetProxy& proxy,
    const IOobjectList& objects
)
{
    const bool syncPar = true;

    // Available fields of type GeoField, sorted order
    const wordList fieldNames(objects.sortedNames<GeoField>());

    // Construct the fields
    PtrList<const GeoField> fields(fieldNames.size());

    label nFields = 0;

    for (const word& fieldName : fieldNames)
    {
        auto tfield =
            getField<GeoField>(proxy, objects, fieldName, syncPar);

        if (tfield)
        {
            fields.set(nFields++, tfield.ptr());
        }
    }

    fields.resize(nFields);
    return fields;
}


// ************************************************************************* //
