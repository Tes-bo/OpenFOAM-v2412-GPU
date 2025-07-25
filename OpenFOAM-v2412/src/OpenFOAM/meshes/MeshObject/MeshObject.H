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

Class
    Foam::MeshObject

Description
    Templated abstract base-class for optional mesh objects used to automate
    their allocation to the mesh database and the mesh-modifier event-loop.

    MeshObject is templated on the type of mesh it is allocated to, the type of
    the mesh object (TopologicalMeshObject, GeometricMeshObject,
    MoveableMeshObject, UpdateableMeshObject) and the type of the actual object
    it is created for.

    Example usage,
    \verbatim
    class leastSquaresVectors
    :
        public MeshObject<fvMesh, MoveableMeshObject, leastSquaresVectors>
    {
    .
    .
    .
        //- Delete the least square vectors when the mesh moves
        virtual bool movePoints();
    };
    \endverbatim

    The MeshObject types:
    - TopologicalMeshObject:
      mesh object to be deleted on topology change
    - GeometricMeshObject:
      mesh object to be deleted on geometry change
    - MoveableMeshObject:
      mesh object to be updated in movePoints
    - UpdateableMeshObject:
      mesh object to be updated in movePoints or updateMesh

Note
    movePoints must be provided for MeshObjects of type MoveableMeshObject
    and both movePoints and updateMesh functions must exist, provided for
    MeshObjects of type UpdateableMeshObject.

SourceFiles
    meshObject.C
    MeshObject.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_MeshObject_H
#define Foam_MeshObject_H

#include "regIOobject.H"
#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                           Class meshObject Declaration
\*---------------------------------------------------------------------------*/

//- The meshObject is a concrete regIOobject to register MeshObject items
class meshObject
:
    public regIOobject
{
public:

    //- Runtime declaration and debug switch
    ClassName("meshObject");


    // Constructors

        //- Construct with given object name on a registry
        meshObject(const word& objName, const objectRegistry& obr);


    // Static Member Functions

        //- Update for mesh motion
        template<class Mesh>
        static void movePoints(objectRegistry& obr);

        //- Update topology using the given map
        template<class Mesh>
        static void updateMesh(objectRegistry& obr, const mapPolyMesh& mpm);

        //- Clear/remove all meshObject of MeshObjectType
        //- via objectRegistry::checkOut()
        template<class Mesh, template<class> class MeshObjectType>
        static void clear(objectRegistry& obr);

        //- Clear all meshObject derived from FromType up to
        //- (but not including) ToType.
        //  Used to clear e.g. all non-updateable meshObjects
        template
        <
            class Mesh,
            template<class> class FromType,
            template<class> class ToType
        >
        static void clearUpto(objectRegistry& obr);
};


/*---------------------------------------------------------------------------*\
                         Class MeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh, template<class> class MeshObjectType, class Type>
class MeshObject
:
    public MeshObjectType<Mesh>
{
protected:

        //- Reference to the mesh
        const Mesh& mesh_;


public:

    // Constructors

        //- Construct with Type::typeName on Mesh
        explicit MeshObject(const Mesh& mesh);

        //- Construct with given object name on Mesh
        MeshObject(const word& objName, const Mesh& mesh);


    //- Destructor
    virtual ~MeshObject() = default;


    // Factory Methods

        //- Get existing or create MeshObject registered with typeName
        template<class... Args>
        static const Type& New(const Mesh& mesh, Args&&... args);

        //- Get existing or create MeshObject with given registration name
        template<class... Args>
        static const Type& New
        (
            const word& objName,
            const Mesh& mesh,
            Args&&... args
        );

        //- Transfer ownership of meshObject to registry.
        static bool Store(std::unique_ptr<Type>&& ptr);

        //- Static destructor using given registration name
        static bool Delete(const word& objName, const Mesh& mesh);

        //- Static destructor using Type::typeName
        static bool Delete(const Mesh& mesh)
        {
            return Delete(Type::typeName, mesh);
        }

        //- Release ownership of meshObject (with given registration name)
        //- from registry. Returns nullptr if not found or not owned.
        static std::unique_ptr<Type> Release
        (
            const word& objName,
            const Mesh& mesh,
            //! optionally perform checkOut() from the registry
            const bool checkout = false
        );


        //- Release ownership of meshObject (with Type::typeName name)
        //- from registry.
        static std::unique_ptr<Type> Release
        (
            const Mesh& mesh,
            //! optionally perform checkOut() from the registry
            const bool checkout = false
        )
        {
            return Release(Type::typeName, mesh, checkout);
        }


    // Member Functions

        //- Reference to the mesh
        const Mesh& mesh() const noexcept
        {
            return mesh_;
        }

        //- Dummy write
        virtual bool writeData(Ostream& os) const
        {
            return true;
        }
};


/*---------------------------------------------------------------------------*\
                    Class TopologicalMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class TopologicalMeshObject
:
    public meshObject
{
public:

    //- Construct from name and instance on registry
    TopologicalMeshObject(const word& objName, const objectRegistry& obr)
    :
        meshObject(objName, obr)
    {}
};


/*---------------------------------------------------------------------------*\
                    Class GeometricMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class GeometricMeshObject
:
    public TopologicalMeshObject<Mesh>
{
public:

    //- Construct from name and instance on registry
    GeometricMeshObject(const word& objName, const objectRegistry& obr)
    :
        TopologicalMeshObject<Mesh>(objName, obr)
    {}
};


/*---------------------------------------------------------------------------*\
                    Class MoveableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class MoveableMeshObject
:
    public GeometricMeshObject<Mesh>
{
public:

    //- Construct from name and instance on registry
    MoveableMeshObject(const word& objName, const objectRegistry& obr)
    :
        GeometricMeshObject<Mesh>(objName, obr)
    {}

    //- Update for mesh motion
    virtual bool movePoints() = 0;
};


/*---------------------------------------------------------------------------*\
                    Class UpdateableMeshObject Declaration
\*---------------------------------------------------------------------------*/

template<class Mesh>
class UpdateableMeshObject
:
    public MoveableMeshObject<Mesh>
{
public:

    //- Construct from name and instance on registry
    UpdateableMeshObject(const word& objName, const objectRegistry& obr)
    :
        MoveableMeshObject<Mesh>(objName, obr)
    {}

    //- Update topology using the given map
    virtual void updateMesh(const mapPolyMesh& mpm) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MeshObject.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
