/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::hexMeshSmootherMotionSolver

Description
    Implementation of hexMeshSmoother (part of extBlockMesh). WIP.

    See https://github.com/Etudes-NG/extBlockMesh
    extBlockMesh Copyright (C) 2014 Etudes-NG

    Quality-based under-relaxation of point smoothing.

Usage
    Example of the motion solver specification in dynamicMeshDict:
    \verbatim
    motionSolver    hexMeshSmoother;
    hexMeshSmootherCoeffs
    {
        //- Number of smoothing iterations
        nPointSmootherIter      10;

        //- Smoother to apply
        pointSmoother           geometricElementTransform;

        //- Any smoother-specific settings
        transformationParameter 0.667;

        //- Underrelax boundary condition
        snapScale               table ((0 0.1) (10 1.0));
    }

SourceFiles
    hexMeshSmootherMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef hexMeshSmootherMotionSolver_H
#define hexMeshSmootherMotionSolver_H

#include "displacementMotionSolver.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

//class searchableSurfaces;
class pointSmoother;

/*---------------------------------------------------------------------------*\
                 Class hexMeshSmootherMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class hexMeshSmootherMotionSolver
:
    public displacementMotionSolver
{
public:

    // Public Data

        //- Enumeration defining the type of attraction
        enum pointType
        {
            INTERIOR = 0,
            SURFACE = 1,
            EDGE = 2,
            POINT = 3
        };


protected:

    // Protected Data

        //- Point smoothing method
        autoPtr<pointSmoother> pointSmoother_;

        //- Number of point smoother iterations per timestep
        const label nPointSmootherIter_;


        // Mesh quality based relaxation of smoothed position

            //- Relaxation factors to use in each iteration
            const scalarList relaxationFactors_;

            //- Per mesh point (including internal) what type
            labelList pointTypes_;

            //- Per mesh point the last used relaxation factor
            labelList relaxationLevel_;

            //- Relaxed point field
            pointField relaxedPoints_;


        ////- Patches that are to be snapped
        //const labelList snapPatches_;
        //
        ////- FaceZones that are to be snapped
        //const labelList snapZones_;

        //- Scaling for snapping
        const autoPtr<Function1<scalar>> snapScale_;

        //- Cached per-point coupled status. For guaranteeing contributions
        //  of coupled points only done once
        const bitSet isMasterPoint_;

        //- Create big primitivePatch for all outside and any features on it
        const autoPtr<indirectPrimitivePatch> bnd0Ptr_;


    // Private Member Functions

        //- Collect all non-constraint (i.e. no processor, cyclic, empty etc)
        //- patches
        static labelList nonConstraintPatches(const polyMesh& mesh);

        //- Create single patch of all supplied faces
        static autoPtr<indirectPrimitivePatch> makePatch
        (
            const polyMesh& mesh,
            const labelList& patchIDs,
            const labelList& zoneIDs,
            const pointField& points0
        );

        //- Check with current points
        void checkMesh
        (
            const pointField& currentPoints,
            const vectorField& fCtrs,
            const vectorField& fAreas,
            const vectorField& cellCtrs,
            const scalarField& cellVols,
            labelHashSet& markedFaces,
            bitSet& markedPoints
        ) const;

        //- Apply current constraints (from pointDisplacement) to supplied
        //  locations
        void constrainDisplacement(pointField& points) const;

        //- Relax the points (supplied in pointDisplacement)
        bool relax
        (
            const scalarList& relaxationFactors,
            const bitSet& pointsToRelax,
            const pointField& initialPoints,
            const pointField& wantedPoints,
            pointField& relaxedPoints,
            labelList& relaxationLevel
        ) const;
        label countPos(const labelList& elems) const;
        labelList countZeroOrPos(const label size, const labelList& lst) const;

        //- Helper: set in bitSet all elements with a certain value
        void select(const labelUList&, const label val, bitSet& isVal) const;

        void laplaceSmooth
        (
            const label type,
            const pointField& initialPoints,
            pointField& newPoints
        ) const;
        void featLaplaceSmooth
        (
            const indirectPrimitivePatch& pp,
            const pointField& initialPoints,
            pointField& newPoints
        ) const;
        // Snap points using boundary evaluation ...
        void snapBoundaryPoints
        (
            const scalar scale,
            const pointField& initialPoints,
            pointField& newPoints
        ) const;

        //- Keep points on empty patches
        void emptyCorrectPoints(pointVectorField& pointDisplacement) const;


public:

    //- Runtime type information
    TypeName("hexMeshSmoother");


    // Constructors

        //- Construct from a polyMesh and an IOdictionary
        hexMeshSmootherMotionSolver
        (
            const polyMesh&,
            const IOdictionary&
        );

        //- Construct from components
        hexMeshSmootherMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict,
            const pointVectorField& pointDisplacement,
            const pointIOField& points0
        );


    //- Destructor
    virtual ~hexMeshSmootherMotionSolver();


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
