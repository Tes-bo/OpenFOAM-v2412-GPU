/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2024 OpenCFD Ltd.
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

#include "MeshObject.H"
#include "objectRegistry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject(const Mesh& mesh)
:
    MeshObjectType<Mesh>(Type::typeName, mesh.thisDb()),
    mesh_(mesh)
{}


template<class Mesh, template<class> class MeshObjectType, class Type>
Foam::MeshObject<Mesh, MeshObjectType, Type>::MeshObject
(
    const word& objName,
    const Mesh& mesh
)
:
    MeshObjectType<Mesh>(objName, mesh.thisDb()),
    mesh_(mesh)
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const Mesh& mesh,
    Args&&... args
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(Type::typeName);

    if (ptr)
    {
        return *ptr;
    }

    if (meshObject::debug)
    {
        Pout<< "MeshObject::New(const " << Mesh::typeName
            << "&, ...) : constructing <" << Type::typeName
            << ">, region=" << mesh.name() << endl;
    }

    ptr = new Type(mesh, std::forward<Args>(args)...);

    regIOobject::store(static_cast<MeshObjectType<Mesh>*>(ptr));

    return *ptr;
}


template<class Mesh, template<class> class MeshObjectType, class Type>
template<class... Args>
const Type& Foam::MeshObject<Mesh, MeshObjectType, Type>::New
(
    const word& objName,
    const Mesh& mesh,
    Args&&... args
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(objName);

    if (ptr)
    {
        return *ptr;
    }

    if (meshObject::debug)
    {
        Pout<< "MeshObject::New('" << objName
            << "', const " << Mesh::typeName
            << "&, ...) : constructing <" << Type::typeName
            << ">, region=" << mesh.name() << endl;
    }

    ptr = new Type(objName, mesh, std::forward<Args>(args)...);

    regIOobject::store(static_cast<MeshObjectType<Mesh>*>(ptr));

    return *ptr;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::MeshObject<Mesh, MeshObjectType, Type>::Delete
(
    const word& objName,
    const Mesh& mesh
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(objName);

    if (ptr)
    {
        if (meshObject::debug)
        {
            Pout<< "MeshObject::Delete() : deleting <" << Type::typeName
                << "> " << objName << endl;
        }

        return mesh.thisDb().checkOut(static_cast<MeshObjectType<Mesh>*>(ptr));
    }

    return false;
}


template<class Mesh, template<class> class MeshObjectType, class Type>
std::unique_ptr<Type> Foam::MeshObject<Mesh, MeshObjectType, Type>::Release
(
    const word& objName,
    const Mesh& mesh,
    const bool checkout
)
{
    Type* ptr =
        mesh.thisDb().objectRegistry::template
        getObjectPtr<Type>(objName);

    std::unique_ptr<Type> released;

    if (ptr)
    {
        auto* casted = static_cast<MeshObjectType<Mesh>*>(ptr);

        if (casted->regIOobject::ownedByRegistry())
        {
            // Release ownership from registry and transfer to unique_ptr
            casted->regIOobject::release();
            released.reset(ptr);

            // Allow removal from the registry (ie, checkOut) but leave its
            // 'registered' status untouched since this is equivalent to
            // IOobject::registerObject().
            //
            // Do not use regIOobject::release(unregister) since this
            // will prevent later re-storing

            if (checkout)
            {
                casted->regIOobject::checkOut();
            }
        }

        if (meshObject::debug)
        {
            Pout<< "MeshObject::Release() : release <" << Type::typeName
                << "> " << objName << ", owned=" << bool(released) << endl;
        }
    }

    return released;
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Mesh, template<class> class MeshObjectType, class Type>
bool Foam::MeshObject<Mesh, MeshObjectType, Type>::Store
(
    std::unique_ptr<Type>&& ptr
)
{
    bool ok = false;

    if (ptr)
    {
        auto* casted = static_cast<MeshObjectType<Mesh>*>(ptr.get());

        ok = casted->regIOobject::store();

        if (ok)
        {
            // Took ownership
            (void) ptr.release();
        }

        if (meshObject::debug)
        {
            Pout<< "MeshObject::Store() : store <" << Type::typeName
                << ">, owned=" << ok << endl;
        }
    }

    return ok;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Mesh>
void Foam::meshObject::movePoints(objectRegistry& obr)
{
    UPtrList<GeometricMeshObject<Mesh>> meshObjects
    (
        obr.sorted<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::movePoints() : moving "
            << meshObjects.size() << " <" << Mesh::typeName
            << "> meshObjects, region=" << obr.name() << endl;
    }

    for (auto& item : meshObjects)
    {
        // isA_constCast<MoveableMeshObject<Mesh>>
        auto* objectPtr = dynamic_cast<MoveableMeshObject<Mesh>*>(&item);

        if (objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Moving " << item.name() << endl;
            }
            objectPtr->movePoints();
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << item.name() << endl;
            }
            obr.checkOut(item);
        }
    }
}


template<class Mesh>
void Foam::meshObject::updateMesh(objectRegistry& obr, const mapPolyMesh& mpm)
{
    UPtrList<GeometricMeshObject<Mesh>> meshObjects
    (
        obr.sorted<GeometricMeshObject<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::updateMesh() : updating "
            << meshObjects.size() << " <" << Mesh::typeName
            << "> meshObjects, region=" << obr.name() << endl;
    }

    for (auto& item : meshObjects)
    {
        // isA_constCast<UpdateableMeshObject<Mesh>>
        auto* objectPtr = dynamic_cast<UpdateableMeshObject<Mesh>*>(&item);

        if (objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Updating " << item.name() << endl;
            }
            objectPtr->updateMesh(mpm);
        }
        else
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << item.name() << endl;
            }
            obr.checkOut(item);
        }
    }
}


template<class Mesh, template<class> class MeshObjectType>
void Foam::meshObject::clear(objectRegistry& obr)
{
    UPtrList<MeshObjectType<Mesh>> meshObjects
    (
        obr.sorted<MeshObjectType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clear() : clearing "
            << meshObjects.size() << " <" << Mesh::typeName
            << "> meshObjects, region=" << obr.name() << endl;
    }

    for (auto& item : meshObjects)
    {
        if (meshObject::debug)
        {
            Pout<< "    Destroying " << item.name() << endl;
        }
        obr.checkOut(item);
    }
}


template
<
    class Mesh,
    template<class> class FromType,
    template<class> class ToType
>
void Foam::meshObject::clearUpto(objectRegistry& obr)
{
    UPtrList<FromType<Mesh>> meshObjects
    (
        obr.sorted<FromType<Mesh>>()
    );

    if (meshObject::debug)
    {
        Pout<< "meshObject::clearUpto() : clearing "
            << meshObjects.size() << " <" << Mesh::typeName
            << "> meshObjects, region=" << obr.name() << endl;
    }

    for (auto& item : meshObjects)
    {
        // isA_constCast<ToType<Mesh>>
        auto* objectPtr = dynamic_cast<ToType<Mesh>*>(&item);

        if (!objectPtr)
        {
            if (meshObject::debug)
            {
                Pout<< "    Destroying " << item.name() << endl;
            }
            obr.checkOut(item);
        }
    }
}


// ************************************************************************* //
