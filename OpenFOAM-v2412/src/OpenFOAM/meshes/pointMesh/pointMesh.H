/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::pointMesh

Description
    Mesh representing a set of points created from polyMesh.

\*---------------------------------------------------------------------------*/

#ifndef Foam_pointMesh_H
#define Foam_pointMesh_H

#include "GeoMesh.H"
#include "MeshObject.H"
#include "polyMesh.H"
#include "pointBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class pointMesh Declaration
\*---------------------------------------------------------------------------*/

class pointMesh
:
    public MeshObject<polyMesh, UpdateableMeshObject, pointMesh>,
    public GeoMesh<polyMesh>
{
    // Private Typedefs

        typedef MeshObject
        <
            polyMesh,
            UpdateableMeshObject,
            pointMesh
        > MeshObject_type;


    // Permanent Data

        //- Boundary mesh
        pointBoundaryMesh boundary_;


    // Private Member Functions

        //- Map all fields
        void mapFields(const mapPolyMesh& mpm);

        //- No copy construct
        pointMesh(const pointMesh&) = delete;

        //- No copy assignment
        void operator=(const pointMesh&) = delete;


public:

    // Public Typedefs

        //- The mesh type
        typedef pointMesh Mesh;

        //- The boundary type associated with the mesh
        typedef pointBoundaryMesh BoundaryMesh;


    // Declare name of the class and its debug switch
    ClassName("pointMesh");

    //- Return the mesh sub-directory name (usually "pointMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from polyMesh
        explicit pointMesh(const polyMesh& pMesh);

        //- Construct from polyMesh and IOobject (used when reading boundary)
        explicit pointMesh(const polyMesh& pMesh, const IOobject& io);

        //- Construct from polyMesh and readOpt. Takes instance, time etc
        //- from polyMesh. Used when reading boundary.
        explicit pointMesh
        (
            const polyMesh& pMesh,
            const IOobjectOption::readOption rOpt
        );


    //- Destructor
    ~pointMesh() = default;


    // Member Functions

        //- Return size. Number of points
        static label size(const Mesh& mesh)
        {
            return mesh.GeoMesh<polyMesh>::mesh_.nPoints();
        }

        //- Return size. Number of points
        label size() const
        {
            return size(*this);
        }

        //- Return reference to boundary mesh
        const pointBoundaryMesh& boundary() const noexcept
        {
            return boundary_;
        }

        //- Is demand-driven parallel info available?
        bool hasGlobalData() const noexcept
        {
            return GeoMesh<polyMesh>::mesh_.hasGlobalData();
        }

        //- Return parallel info (demand-driven)
        const globalMeshData& globalData() const
        {
            return GeoMesh<polyMesh>::mesh_.globalData();
        }

        //- Return database. For now is its polyMesh.
        const objectRegistry& thisDb() const
        {
            return GeoMesh<polyMesh>::mesh_.thisDb();
        }

        //- Return Time from polyMesh.
        const Time& time() const
        {
            return GeoMesh<polyMesh>::mesh_.time();
        }

        //- Set the instance for mesh files
        void setInstance
        (
            const fileName& instance,
            const IOobjectOption::writeOption wOpt = IOobject::AUTO_WRITE
        );


    // Volume Mesh

        //- Return boundary mesh for underlying volume mesh
        const polyBoundaryMesh& boundaryMesh() const
        {
            return GeoMesh<polyMesh>::mesh_.boundaryMesh();
        }


    // Mesh Motion

        //- Move points
        bool movePoints();

        //- Update the mesh corresponding to given map
        void updateMesh(const mapPolyMesh& mpm);


    // Member Operators

        bool operator!=(const pointMesh& pm) const
        {
            return &pm != this;
        }

        bool operator==(const pointMesh& pm) const
        {
            return &pm == this;
        }


    // Write

        //- Write
        virtual bool writeObject
        (
            IOstreamOption streamOpt,
            const bool writeOnProc = true
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
