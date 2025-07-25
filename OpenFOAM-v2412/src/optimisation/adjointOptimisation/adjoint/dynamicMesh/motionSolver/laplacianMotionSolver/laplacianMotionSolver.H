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
    Foam::laplacianMotionSolver

Description
    Similar to velocityLaplacian but iteratively solves the mesh displacement
    PDEs to account for non-orthogonality.
    The boundary displacement is set as a boundary condition
    on pointMotionU; the latter is generated automatically if not found.
    Assumes uniform diffusivity

SourceFiles
    laplacianMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef laplacianMotionSolver_H
#define laplacianMotionSolver_H

#include "velocityMotionSolver.H"
#include "fvMotionSolver.H"
#include "volPointInterpolation.H"
#include "polyMesh.H"
#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
class motionInterpolation;
class motionDiffusivity;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                   Class laplacianMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class laplacianMotionSolver
:
    public motionSolver,
    public fvMotionSolver
{
protected:

    // Protected data

        mutable pointVectorField pointMotionU_;
        volVectorField cellMotionU_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Diffusivity used to control the motion
        autoPtr<motionDiffusivity> diffusivityPtr_;

        //- Number of laplacian iterations per solution step
        label nIters_;

        //- Residual threshold
        scalar tolerance_;


private:


    // Private Member Functions

        //- No copy construct
        laplacianMotionSolver(const laplacianMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const laplacianMotionSolver&) = delete;

public:

    //- Runtime type information
    TypeName("laplacianMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        laplacianMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~laplacianMotionSolver() = default;


    // Member Functions

        //- Get const and non-const references to pointMotionU
        inline pointVectorField& pointMotionU();
        inline const pointVectorField& pointMotionU() const;

        //- Get const and non-const references to cellMotionU
        inline volVectorField& cellMotionU();
        inline const volVectorField& cellMotionU() const;

        //- Return point location obtained from the current motion field
        virtual tmp<pointField> curPoints() const;

        //- Solve for motion
        virtual void solve();

        //- Set boundary conditions of cellMotionU based on pointMotionU.
        //  Avoiding the use of the cellMotionFvPatchField bc
        //  due to the use of the registry in order to grab pointMotionU, which
        //  may give problems if multiple objects of this class are constructed
        //  at the same time
        void setBoundaryConditions();

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "laplacianMotionSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
