/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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
    Foam::faMeshesRegistry

Description
    A MeshObject registry on top of a polyMesh that provides
    an objectRegistry for holding multiple faMesh objects.

    Directory structure for fields/meshes (default region):
    \verbatim
    .
    |-- constant
    |   `-- finite-area
    |       `-- faMesh
    |           |-- faceLabels
    |           |-- ...
    |           `-- faBoundary
    |-- system
    |   `-- finite-area
    |       |-- faMeshDefiniton
    |       |-- faSchemes
    |       `-- faSolution
    `-- instance
        `-- finite-area
            |-- U
            |-- p
            `-- ...
    \endverbatim

    Directory structure for fields/meshes (multi-regions):
    \verbatim
    .
    |-- constant
    |   `-- finite-area
    |       |-- regionName1
    |       |   `-- faMesh
    |       |       `-- ...
    |       `-- regionNameN
    |           `-- faMesh
    |               `-- ...
    |-- system
    |   `-- finite-area
    |       |-- regionName1
    |       |   |-- faMeshDefiniton
    |       |   |-- faSchemes
    |       |   `-- faSolution
    |       `-- regionNameN
    |           `-- ...
    |
    `-- instance
        `-- finite-area
            |-- regionName1
            |   `-- ...
            `-- regionNameN
                `-- ...
    \endverbatim

SourceFiles
    faMeshesRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_faMeshesRegistry_H
#define Foam_faMeshesRegistry_H

#include "MeshObject.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class faMesh;

/*---------------------------------------------------------------------------*\
                      Class faMeshesRegistry Declaration
\*---------------------------------------------------------------------------*/

class faMeshesRegistry
:
    public MeshObject<polyMesh, UpdateableMeshObject, faMeshesRegistry>
{
    // Data Types

        //- The MeshObject type
        typedef MeshObject
        <
            polyMesh,
            UpdateableMeshObject,
            faMeshesRegistry
        > MeshObject_type;


    // Private Member Data

        //- The sub-registry of finite-area objects (meshes),
        //- anchored on the polyMesh parent. Name: "finite-area",
        objectRegistry objects_;


public:

    //- Runtime type information.
    TypeName("faMeshesRegistry");

    // Generated Methods

        //- No copy construct
        faMeshesRegistry(const faMeshesRegistry&) = delete;

        //- No copy assignment
        void operator=(const faMeshesRegistry&) = delete;


    // Constructors

        //- Construct as singleton on the polyMesh registry
        explicit faMeshesRegistry(const polyMesh& mesh);


    // Database

        //- Return the object registry
        const objectRegistry& thisDb() const noexcept
        {
            return objects_;
        }

        //- The polyMesh reference
        const polyMesh& mesh() const noexcept
        {
            return MeshObject_type::mesh();
        }


    // Topological Change

        //- Is mesh moving - ie, is polyMesh moving
        bool moving() const { return MeshObject_type::mesh().moving(); }

        //- Update after mesh motion
        virtual bool movePoints();

        //- Update after topo change
        virtual void updateMesh(const mapPolyMesh& mpm);


    // Write

        //- Write items (eg, faMesh) held in the registry
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
