/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::elasticityMotionSolver

Description
    Mesh deformation based on the linear elasticity equations.
    The boundary displacement is set as a boundary condition
    on the pointMotionU field.

    Reference:
    \verbatim
        Dwight, R. P., (2009).
        Robust Mesh Deformation using the Linear Elasticity Equations.
        Computational Fluid Dynamics 2006, 401-406.
        Springer Berlin Heidelberg.
    \endverbatim

SourceFiles
    elasticityMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef elasticityMotionSolver_H
#define elasticityMotionSolver_H

#include "motionSolver.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class elasticityMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class elasticityMotionSolver
:
    public motionSolver
{
protected:

    // Protected data

        //- Since the mesh deformation is broken down to multiple steps,
        //- mesh points need to be moved here.
        //- Hence, the non-const mesh reference
        fvMesh& fvMesh_;
        pointVectorField pointMotionU_;
        volVectorField cellMotionU_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Intermediate steps to solve the PDEs
        label nSteps_;

        //- Number of laplacian iterations per solution step
        label nIters_;

        //- Residual threshold
        scalar tolerance_;


    // Protected Member Functions

        //- Set boundary conditions of cellMotionU based on pointMotionU.
        //  Avoiding the use of the cellMotionFvPatchField bc
        //  due to the use of the registry in order to grab pointMotionU, which
        //  may give problems if multiple objects of this class are constructed
        //  at the same time
        void setBoundaryConditions();


private:

    // Private Member Functions

        //- No copy construct
        elasticityMotionSolver(const elasticityMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const elasticityMotionSolver&) = delete;


public:

    //- Runtime type information
    TypeName("elasticityMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        elasticityMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~elasticityMotionSolver() = default;


    // Member Functions

        //- Get const and non-const references to pointMotionU
        inline pointVectorField& pointMotionU();
        inline const pointVectorField& pointMotionU() const;

        //- Get const and non-const references to cellMotionU
        inline volVectorField& cellMotionU();
        inline const volVectorField& cellMotionU() const;

        //- Return point location. Mesh is actually moved in solve()
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion.
        //  Does the actual mesh displacement here, since it is broken down
        //  into multiple steps
        virtual void solve();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "elasticityMotionSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
