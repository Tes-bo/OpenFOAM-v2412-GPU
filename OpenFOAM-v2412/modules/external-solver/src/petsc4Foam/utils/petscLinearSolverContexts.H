/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::petscLinearSolverContexts

Description
    A mesh-based collection linear solver contexts for petsc

SourceFiles
    petscLinearSolverContexts.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_module_petscLinearSolverContexts_H
#define Foam_module_petscLinearSolverContexts_H

#include "fvMesh.H"
#include "MeshObject.H"
#include "HashPtrTable.H"
#include "petscLinearSolverContext.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class petscLinearSolverContexts Declaration
\*---------------------------------------------------------------------------*/

class petscLinearSolverContexts
:
    public MeshObject
    <
        fvMesh,
        TopologicalMeshObject,
        petscLinearSolverContexts
    >
{
    // Private Data

        //- Linear solver contexts
        mutable HashPtrTable<petscLinearSolverContext> contexts_;


public:

    //- Runtime information
    TypeName("petscLinearSolverContext");


    // Constructors

        //- Construct from mesh for given field
        explicit petscLinearSolverContexts(const fvMesh& mesh)
        :
            MeshObject
            <
                fvMesh,
                Foam::TopologicalMeshObject,
                petscLinearSolverContexts
            >(mesh),
            contexts_()
        {}


    //- Destructor
    virtual ~petscLinearSolverContexts() = default;


    // Member Functions

        //- Linear solver contexts (mutable)
        HashPtrTable<petscLinearSolverContext>& contexts() const
        {
            return contexts_;
        }


        //- Test if named context exists
        bool haveContext(const word& name) const
        {
            return contexts_.found(name);
        }


        //- Get existing or a new context by name
        petscLinearSolverContext& getContext(const word& name) const
        {
            auto iter = contexts().find(name);

            const bool cached = iter.found();

            if (!iter.found())
            {
                contexts().set(name, new petscLinearSolverContext{});
                iter = contexts().find(name);
            }

            DebugInfo
                << "cache-" << (cached ? "hit" : "miss")
                << " PETSc linear solver context "
                << name << " mesh: " << mesh().name() << nl;

            return *(*iter);
        }


        //- Remove context by name
        bool removeContext(const word& name) const
        {
            return contexts().erase(name);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //



