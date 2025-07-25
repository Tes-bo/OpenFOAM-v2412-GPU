/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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
    Foam::pLaplacianMotionSolver

Description
    Similar to velocityLaplacian but with a variable diffusivity, based
    on the gradient of the displacement.
    The boundary displacement is set as a boundary condition
    on pointMotionU; the latter is generated automatically if not found.

SourceFiles
    pLaplacianMotionSolver.C

\*---------------------------------------------------------------------------*/

#ifndef pLaplacianMotionSolver_H
#define pLaplacianMotionSolver_H

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
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                    Class pLaplacianMotionSolver Declaration
\*---------------------------------------------------------------------------*/

class pLaplacianMotionSolver
:
    public motionSolver,
    public fvMotionSolver
{
protected:

    // Protected data

        // Use a fixedValue boundary condition for the moving patches?  When
        // using this motionSolver in an optimisation context, we usually go
        // from point movement to face movement, solve the grid dispalcement
        // PDE and then interpolate the movement back to the points. If the
        // boundary conditions of pointMotionU are fixedValue, the initial
        // values will be retained, otherwise they will be overwritten by the
        // face-to-point interpolation. The latter is usually beneficial for
        // the resulting mesh quality but does not give us the exact geometry.
        bool useFixedValuePointMotionUBCs_;

        mutable pointVectorField pointMotionU_;
        volVectorField cellMotionU_;

        //- Interpolation used to transfer cell displacement to the points
        autoPtr<motionInterpolation> interpolationPtr_;

        //- Number of pLaplacian iterations per solution step
        label nIters_;

        //- Residual threshold
        scalar tolerance_;

        //- Residual threshold for intermediate exponents
        scalar toleranceIntermediate_;

        //- Exponent defining the order or the p-Laplacian
        label exponent_;


private:


    // Private Member Functions

        //- No copy construct
        pLaplacianMotionSolver(const pLaplacianMotionSolver&) = delete;

        //- No copy assignment
        void operator=(const pLaplacianMotionSolver&) = delete;

public:

    //- Runtime type information
    TypeName("pLaplacianMotionSolver");


    // Constructors

        //- Construct from mesh and dictionary
        pLaplacianMotionSolver
        (
            const polyMesh& mesh,
            const IOdictionary& dict
        );


    //- Destructor
    virtual ~pLaplacianMotionSolver() = default;


    // Member Functions

        //- Set boundary conditions of cellMotionU based on pointMotionU.
        //  Avoiding the use of the cellMotionFvPatchField bc
        //  due to the use of the registry in order to grab pointMotionU, which
        //  may give problems if multiple objects of this class are constructed
        //  at the same time
        void setBoundaryConditions();

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

        //- Update local data for geometry changes
        virtual void movePoints(const pointField&);

        //- Update the mesh corresponding to given map
        virtual void updateMesh(const mapPolyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pLaplacianMotionSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
