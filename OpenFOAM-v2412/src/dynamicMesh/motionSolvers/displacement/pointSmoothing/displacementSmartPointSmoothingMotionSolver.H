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
    Foam::displacementSmartPointSmoothingMotionSolver

Description
    Quality-based under-relaxation for run-time selectable point smoothing. WIP.

Usage
    \table
        Property     | Description                | Required | Default value
        untangler    | pointSmoother for untangling | yes       | 
        untangleQ    | quality below which untangling is applied | yes  | 

        pointSmoother | pointSmoother in 'normal' mode | yes
        minQ          | quality below which pointSmoother is applied | yes
        nPointSmootherIter | max number of iterations
    \endtable

    Example of the motion solver specification in dynamicMeshDict:
    \verbatim
    motionSolver    displacementSmartPointSmoothing;
    displacementSmartPointSmoothingCoeffs
    {
        //- Overall max number of smoothing iterations
        nPointSmootherIter      10;

        //- If any faces have quality below untangleQ apply untangler
        untangleQ               0.001;
        untangler               laplacian;

        //- If any faces have quality below minQ apply 'normal' smoother
        minQ                    0.9;
        pointSmoother           geometricElementTransform;
        transformationParameter 0.667;
    }
    \endverbatim

SourceFiles
    displacementSmartPointSmoothingMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef displacementSmartPointSmoothingMotionSolver_H
#define displacementSmartPointSmoothingMotionSolver_H

#include "displacementMotionSolver.H"
#include "pointSmoother.H"
#include "polyMeshGeometry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
         Class displacementSmartPointSmoothingMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class displacementSmartPointSmoothingMotionSolver
:
    public displacementMotionSolver
{
protected:

    // Protected Data

        //- Part-updatable mesh geometry
        polyMeshGeometry meshGeometry_;

        //- Point untangler method
        autoPtr<pointSmoother> pointUntangler_;

        //- Minimum allowed quality
        const scalar untangleQ_;

        //- Minimum allowed quality
        const scalar minQ_;

        //- Point smoothing method
        autoPtr<pointSmoother> pointSmoother_;

        //- Number of point smoother iterations per timestep
        const label nPointSmootherIter_;


        // Mesh quality based relaxation of smoothed position

            //- Relaxation factors to use in each iteration
            scalarList relaxationFactors_;

            //- Relaxed point field
            pointField relaxedPoints_;

            //- Set of the faces which are to be moved
            labelHashSet facesToMove_;

            //- Mesh quality dictionary
            dictionary meshQualityDict_;


    // Private Member Functions

        //- Mark affected faces
        void markAffectedFaces
        (
            const labelHashSet& changedFaces,
            labelHashSet& affectedFaces
        );

        //- Relax the points
        bool relax();

        //- Set all the faces to be moved
        void virtual setFacesToMove(const dictionary&);

        //- Handle 2D & empty bcs. Assume in both cases the starting mesh
        //      - has all edges aligned with 3rd dimension
        //      - is on planes of the empty patches
        void emptyCorrectPoints(pointVectorField& pointDisplacement);


public:

    //- Runtime type information
    TypeName("displacementSmartPointSmoothing");


    // Constructors

        //- Construct from a polyMesh and an IOdictionary
        displacementSmartPointSmoothingMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        displacementSmartPointSmoothingMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    virtual ~displacementSmartPointSmoothingMotionSolver() = default;


    // Member Functions

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
