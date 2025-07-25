/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
template<class... Args>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New_impl
(
    IOobjectOption::registerOption regOpt,
    const word& name,
    const Mesh& mesh,
    Args&&... args
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        // == mesh.thisDb().newIOobject(name)
        IOobject
        (
            name,
            mesh.thisDb().time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        mesh,
        std::forward<Args>(args)...
    );

    // Registration
    {
        if (IOobjectOption::REGISTER == regOpt)
        {
            ptr->checkIn();
        }
        else if
        (
            // LEGACY_REGISTER: detect if caching is desired
            (IOobjectOption::LEGACY_REGISTER == regOpt)
         && ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        iField,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        iField,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        std::move(iField),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        std::move(iField),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField,
    const PtrList<PatchField<Type>>& pflds
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        iField,
        pflds
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField,
    const PtrList<PatchField<Type>>& pflds
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        iField,
        pflds
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField,
    const PtrList<PatchField<Type>>& pflds
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        std::move(iField),
        pflds
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField,
    const PtrList<PatchField<Type>>& pflds
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        std::move(iField),
        pflds
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        value,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        value,
        dims,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        value,
        dims,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        value,
        dims,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        regOpt,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        regOpt,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions(),
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            tfld().instance(),
            tfld().local(),
            tfld().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        tfld
    );

    // Registration
    {
        if (IOobjectOption::REGISTER == regOpt)
        {
            ptr->checkIn();
        }
        else if
        (
            // LEGACY_REGISTER: detect if caching is desired
            (IOobjectOption::LEGACY_REGISTER == regOpt)
         && ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        IOobjectOption::LEGACY_REGISTER,
        tfld
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld,
    const word& patchFieldType
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            tfld().instance(),
            tfld().local(),
            tfld().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        tfld,
        patchFieldType
    );

    // Registration
    {
        if (IOobjectOption::REGISTER == regOpt)
        {
            ptr->checkIn();
        }
        else if
        (
            // LEGACY_REGISTER: detect if caching is desired
            (IOobjectOption::LEGACY_REGISTER == regOpt)
         && ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld,
    const word& patchFieldType
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        IOobjectOption::LEGACY_REGISTER,
        tfld,
        patchFieldType
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            tfld().instance(),
            tfld().local(),
            tfld().db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        tfld,
        patchFieldTypes,
        actualPatchTypes
    );

    // Registration
    {
        if (IOobjectOption::REGISTER == regOpt)
        {
            ptr->checkIn();
        }
        else if
        (
            // LEGACY_REGISTER: detect if caching is desired
            (IOobjectOption::LEGACY_REGISTER == regOpt)
         && ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


template<class Type, template<class> class PatchField, class GeoMesh>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const word& name,
    const tmp<GeometricField<Type, PatchField, GeoMesh>>& tfld,
    const wordList& patchFieldTypes,
    const wordList& actualPatchTypes
)
{
    return GeometricField<Type, PatchField, GeoMesh>::New
    (
        name,
        IOobjectOption::LEGACY_REGISTER,
        tfld,
        patchFieldTypes,
        actualPatchTypes
    );
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const GeometricField<AnyType, PatchField, GeoMesh>& fld,
    const word& name,
    const dimensionSet& dims,
    const word& patchFieldType
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            fld.instance(),
            fld.db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        fld.mesh(),
        dims,
        patchFieldType
    );

    // Registration
    {
        if
        (
            // LEGACY_REGISTER: detect if caching is desired
            ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


template<class Type, template<class> class PatchField, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::GeometricField<Type, PatchField, GeoMesh>>
Foam::GeometricField<Type, PatchField, GeoMesh>::New
(
    const GeometricField<AnyType, PatchField, GeoMesh>& fld,
    const word& name,
    const dimensioned<Type>& dt,
    const word& patchFieldType
)
{
    auto ptr = tmp<GeometricField<Type, PatchField, GeoMesh>>::New
    (
        IOobject
        (
            name,
            fld.instance(),
            fld.db(),
            IOobjectOption::NO_READ,
            IOobjectOption::NO_WRITE,
            IOobjectOption::NO_REGISTER
        ),
        fld.mesh(),
        dt.value(),
        dt.dimensions(),
        patchFieldType
    );

    // Registration
    {
        if
        (
            // LEGACY_REGISTER: detect if caching is desired
            ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


// ************************************************************************* //
