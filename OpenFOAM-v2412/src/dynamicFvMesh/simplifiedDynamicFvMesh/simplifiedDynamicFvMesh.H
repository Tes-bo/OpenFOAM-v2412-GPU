/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::simplifiedDynamicFvMesh

Description
    Functions to generate simplified finite volume meshes

SourceFiles
    simplifiedDynamicFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef simplifiedDynamicFvMesh_H
#define simplifiedDynamicFvMesh_H

#include "columnFvMesh.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dynamicFvMesh;

namespace simplifiedMeshes
{

/*---------------------------------------------------------------------------*\
                   Class simplifiedDynamicFvMesh Declaration
\*---------------------------------------------------------------------------*/

class simplifiedDynamicFvMeshBase
{
public:

    //- Runtime type information
    TypeName("simplifiedDynamicFvMeshBase");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            dynamicFvMesh,
            time,
            (
                const Time& runTime,
                const word& regionName
            ),
            (runTime, regionName)
        );


    // Selectors

        //- Return a reference to the selected simplified mesh
        static autoPtr<dynamicFvMesh> New(const IOobject& io);


    //- Constructor
    simplifiedDynamicFvMeshBase()
    {}

    //- Destructor
    virtual ~simplifiedDynamicFvMeshBase() = default;
};


template<class DynamicMeshType>
class SimplifiedDynamicFvMesh
:
    public simplifiedDynamicFvMeshBase,
    public columnFvMeshInfo,
    public DynamicMeshType
{

public:

    ClassNameNoDebug(DynamicMeshType::typeName_.c_str());

    //- Constructor
    SimplifiedDynamicFvMesh(const Time& runTime, const word& regionName);

    //- Update the mesh for both mesh motion and topology change
    virtual bool update()
    {
        // No updates performed
        return false;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace simplifiedMeshes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "simplifiedDynamicFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define createProxyDynamicFvMesh(Type)                                        \
                                                                              \
typedef simplifiedMeshes::SimplifiedDynamicFvMesh<Type> simplified##Type;     \
                                                                              \
template<>                                                                    \
const word simplified##Type::typeName(Type::typeName_());                     \
                                                                              \
namespace simplifiedMeshes                                                    \
{                                                                             \
    addToRunTimeSelectionTable                                                \
    (                                                                         \
        simplifiedDynamicFvMeshBase,                                          \
        simplified##Type,                                                     \
        time                                                                  \
    );                                                                        \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
