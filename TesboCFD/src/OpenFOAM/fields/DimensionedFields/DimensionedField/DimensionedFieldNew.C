/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

template<class Type, class GeoMesh>
template<class... Args>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New_impl
(
    IOobjectOption::registerOption regOpt,
    const word& name,
    const Mesh& mesh,
    Args&&... args
)
{
    auto ptr = tmp<DimensionedField<Type, GeoMesh>>::New
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

template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        iField
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    const Field<Type>& iField
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        iField
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        std::move(iField)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims,
    Field<Type>&& iField
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        std::move(iField)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensionSet& dims
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        dims,
        false  // No checkIOFlags (is NO_READ anyhow)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensionSet& dims
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        dims,
        false  // No checkIOFlags (is NO_READ anyhow)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        regOpt,
        name,
        mesh,
        value,
        dims,
        false  // No checkIOFlags (is NO_READ anyhow)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const Type& value,
    const dimensionSet& dims
)
{
    return DimensionedField<Type, GeoMesh>::New_impl
    (
        IOobjectOption::LEGACY_REGISTER,
        name,
        mesh,
        value,
        dims,
        false  // No checkIOFlags (is NO_READ anyhow)
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    return DimensionedField<Type, GeoMesh>::New
    (
        name,
        regOpt,
        mesh,
        dt.value(),
        dt.dimensions()
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const Mesh& mesh,
    const dimensioned<Type>& dt
)
{
    return DimensionedField<Type, GeoMesh>::New
    (
        name,
        mesh,
        dt.value(),
        dt.dimensions()
    );
}


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    IOobjectOption::registerOption regOpt,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    auto ptr = tmp<DimensionedField<Type, GeoMesh>>::New
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


template<class Type, class GeoMesh>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const word& name,
    const tmp<DimensionedField<Type, GeoMesh>>& tfld
)
{
    return DimensionedField<Type, GeoMesh>::New
    (
        name,
        IOobjectOption::LEGACY_REGISTER,
        tfld
    );
}


template<class Type, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const DimensionedField<AnyType, GeoMesh>& fld,
    const word& name,
    const dimensionSet& dims
)
{
    auto ptr = tmp<DimensionedField<Type, GeoMesh>>::New
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
        dims
    );

    // Registration
    {
        if
        (
            ptr->db().is_cacheTemporaryObject(ptr.get())
        )
        {
            ptr.protect(true);
            ptr->checkIn();
        }
    }
    return ptr;
}


template<class Type, class GeoMesh>
template<class AnyType>
Foam::tmp<Foam::DimensionedField<Type, GeoMesh>>
Foam::DimensionedField<Type, GeoMesh>::New
(
    const DimensionedField<AnyType, GeoMesh>& fld,
    const word& name,
    const dimensioned<Type>& dt
)
{
    auto ptr = tmp<DimensionedField<Type, GeoMesh>>::New
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
        dt.dimensions()
    );

    // Registration
    {
        if
        (
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
