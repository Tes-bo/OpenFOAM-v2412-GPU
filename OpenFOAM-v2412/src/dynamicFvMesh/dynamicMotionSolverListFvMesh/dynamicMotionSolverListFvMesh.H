/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
    Copyright (C) 2020,2022 OpenCFD Ltd.
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
    Foam::dynamicMotionSolverListFvMesh

Description
    Dynamic mesh able to handle multiple motion solvers.
    NOTE: If the word entry "solvers" is not found it falls back to a single
    motion solver behavior.

SourceFiles
    dynamicMotionSolverListFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverListFvMesh_H
#define dynamicMotionSolverListFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                Class dynamicMotionSolverListFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverListFvMesh
:
    public dynamicFvMesh
{
    // Private data

        PtrList<motionSolver> motionSolvers_;


    // Private Member Functions

        //- No copy construct
        dynamicMotionSolverListFvMesh
        (
            const dynamicMotionSolverListFvMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const dynamicMotionSolverListFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverListFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverListFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    ~dynamicMotionSolverListFvMesh();


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Initialise all non-demand-driven data. Used for initialising
        //  without default motionSolver (used by dynamicRefineFvMesh)
        virtual bool init(const bool doInit, const bool mandatory);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();

        //- Map all fields in time using given map. Triggered by topo change
        virtual void mapFields(const mapPolyMesh& mpm);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
