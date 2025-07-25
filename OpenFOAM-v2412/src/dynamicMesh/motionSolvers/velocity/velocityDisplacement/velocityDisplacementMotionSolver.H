/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::velocityDisplacementMotionSolver

Description
    Mesh motion solver for a polyMesh. Wraps a displacement motion solver in a
    velocity motion solver.

SourceFiles
    velocityDisplacementMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef velocityDisplacementMotionSolver_H
#define velocityDisplacementMotionSolver_H

#include "velocityMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class displacementMotionSolver;

/*---------------------------------------------------------------------------*\
              Class velocityDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class velocityDisplacementMotionSolver
:
    public velocityMotionSolver
{
    // Private data

        //- Displacement motion solver
        autoPtr<displacementMotionSolver> displacementMotionSolverPtr_;


    // Private Member Functions

        //- No copy construct
        velocityDisplacementMotionSolver
        (
            const velocityDisplacementMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const velocityDisplacementMotionSolver&) = delete;

        //- Get the boundary condition types for the point displacement
        wordList pointDisplacementBoundaryTypes() const;


public:

    //- Runtime type information
    TypeName("velocityDisplacement");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        velocityDisplacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );


    //- Destructor
    ~velocityDisplacementMotionSolver();


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update geometry
        virtual void movePoints(const pointField&);

        //- Update topology
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
