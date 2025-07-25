/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2020 OpenCFD Ltd
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
    Foam::dynamicMotionSolverTopoFvMesh

Description
    Sample topoChangerFvMesh that moves an object using a motion solver

SourceFiles
    dynamicMotionSolverTopoFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverTopoFvMesh_H
#define dynamicMotionSolverTopoFvMesh_H

#include "topoChangerFvMesh.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                      Class dynamicMotionSolverTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary motionDict_;

        //- Motion solver
        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- No copy construct
        dynamicMotionSolverTopoFvMesh
        (
            const dynamicMotionSolverTopoFvMesh&
        ) = delete;

        //- No copy assignment
        void operator=(const dynamicMotionSolverTopoFvMesh&) = delete;


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit dynamicMotionSolverTopoFvMesh
        (
            const IOobject& io,
            const bool doInit=true
        );


    //- Destructor
    virtual ~dynamicMotionSolverTopoFvMesh();


    // Member Functions

        //- Initialise all non-demand-driven data
        virtual bool init(const bool doInit);

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
