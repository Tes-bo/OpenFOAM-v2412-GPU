/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2024 OpenCFD Ltd.
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
    Foam::multiDisplacementMotionSolver

Group
    grpMeshMotionSolvers

Description
    Mesh motion solver for a polyMesh. Applies multiple (displacement) motion
    solvers in order.

    Not very efficient : all displacementMotionSolvers store a copy
    of the initial points (points0) and the displacement (pointDisplacement
    or also cellDisplacement).
    Used to combine large-scale, implicit displacement smoothing (e.g.
    displacementLaplacian) with point smoothing.

Usage
    Example of the dynamicMeshDict specification:
    \verbatim
    motionSolver    multiDisplacement;
    solvers
    {
        // Solve finite volume laplacian to efficiently smooth displacement
        // (not point locations)
        displacementLaplacian
        {
            motionSolver        displacementLaplacian;
            diffusivity         uniform;
        }

        // Apply few iterations of smoothing of point locations
        displacementPointSmoothing
        {
            motionSolver        displacementPointSmoothing;
            pointSmoother       laplacian;
            nPointSmootherIter  10;
        }
    }
    \endverbatim

Note
    When using displacementLaplacian: the default behaviour for the
    cellDisplacement is to apply fixed value boundary conditions (by averaging
    point values) only to those pointDisplacement boundary conditions that
    are fixed value. Quite a few point boundary conditions (e.g. surfaceSlip,
    edgeSlip) are not so require an explicitly provided cellDisplacement
    field with 'cellMotion' boundary conditions for those patches.

SourceFiles
    isplacementMultiMotionSolver.C
\*----------------------------------------------------------------------------*/

#ifndef Foam_multiDisplacementMotionSolver_H
#define Foam_multiDisplacementMotionSolver_H

#include "displacementMotionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class multiDisplacementMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class multiDisplacementMotionSolver
:
    public displacementMotionSolver
{
    // Private data

        //- Current points
        pointField curPoints_;

        //- List of motion solvers
        PtrList<displacementMotionSolver> motionSolvers_;


    // Private Member Functions

        //- No copy construct
        multiDisplacementMotionSolver
        (
            const multiDisplacementMotionSolver&
        ) = delete;

        //- No copy assignment
        void operator=(const multiDisplacementMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("multiDisplacement");


    // Constructors

        //- Construct from polyMesh and IOdictionary
        multiDisplacementMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        multiDisplacementMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    ~multiDisplacementMotionSolver() = default;


    // Member Functions

        //- Provide current points for motion
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //-  Update local data for topology changes
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
