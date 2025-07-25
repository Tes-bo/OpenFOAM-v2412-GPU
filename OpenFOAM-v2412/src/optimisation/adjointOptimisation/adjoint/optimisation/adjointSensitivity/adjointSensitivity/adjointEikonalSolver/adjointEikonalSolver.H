/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
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
    Foam::adjointEikonalSolver

Description
    Solver of the adjoint to the eikonal PDE

    Reference:
    \verbatim
        For the development of the adjoint eikonal PDE and its boundary
        conditions

            Papoutsis-Kiachagias, E. M., & Giannakoglou, K. C. (2014).
            Continuous Adjoint Methods for Turbulent Flows, Applied to Shape
            and Topology Optimization: Industrial Applications.
            Archives of Computational Methods in Engineering, 23(2), 255-299.
            http://doi.org/10.1007/s11831-014-9141-9
    \endverbatim

    To be as consistent as possible, it is recommended to use the
    advectionDiffusion wallDist method in fvSchemes, instead of the more widely
    used meshWave

    Example of the adjointEikonalSolver specification in optimisationDict:
    \verbatim
        optimisation
        {
            sensitivities
            {
                includeDistance true;
                adjointEikonalSolver
                {
                    // epsilon should be the same as the one used
                    // in fvSchemes/wallDist/advectionDiffusionCoeffs
                    epsilon   0.1;
                    iters     1000;
                    tolerance 1e-6;
                }
            }
        }
    \endverbatim

    Example of the entries in fvSchemes:
    \verbatim
        divSchemes
        {
            .
            .
            // avoid bounded schemes since yPhi is not conservative
            div(-yPhi,da) Gauss linearUpwind grad(da);
            .
            .
        }
        laplacianSchemes
        {
            .
            .
            laplacian(yWall,da) Gauss linear corrected;
            .
            .
        }
    \endverbatim

    Also, the solver specification and a relaxation factor for da are required
    in fvSolution

    \verbatim
        da
        {
            solver          PBiCGStab;
            preconditioner  DILU;
            tolerance       1e-9;
            relTol          0.1;
        }

        relaxationFactors
        {
            equations
            {
                .
                .
                da           0.5;
                .
                .
            }
        }
    \endverbatim

See also
    Foam::patchDistMethod::advectionDiffusion
    Foam::wallDist


SourceFiles
    adjointEikonalSolver.C

\*---------------------------------------------------------------------------*/

#ifndef adjointEikonalSolver_H
#define adjointEikonalSolver_H

#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "fvMesh.H"
#include "calculatedFvPatchField.H"
#include "createZeroField.H"
#include "boundaryFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declaration
class adjointSolver;

/*---------------------------------------------------------------------------*\
                     Class adjointEikonalSolver Declaration
\*---------------------------------------------------------------------------*/

class adjointEikonalSolver
{
private:

    // Private Member Functions

        //- No copy construct
        adjointEikonalSolver(const adjointEikonalSolver&) = delete;

        //- No copy assignment
        void operator=(const adjointEikonalSolver&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;

        dictionary dict_;

        adjointSolver& adjointSolver_;

        const labelHashSet& sensitivityPatchIDs_;

        label nEikonalIters_;

        scalar tolerance_;

        scalar epsilon_;

        labelHashSet wallPatchIDs_;

        volScalarField da_;

        volScalarField source_;

        //- Wall face sens w.r.t. (x,y.z)
        autoPtr<boundaryVectorField> distanceSensPtr_;


    // Protected functions

        //- Return the boundary condition types for da
        wordList patchTypes() const;

        //- Compute convecting velocity
        tmp<surfaceScalarField> computeYPhi();

        //- Read options each time a new solution is found
        void read();



public:

    //- Runtime type information
    TypeName("adjointEikonalSolver");


    // Constructors

        //- Construct from components
        adjointEikonalSolver
        (
            const fvMesh& mesh,
            const dictionary& dict,
            adjointSolver& adjointSolver,
            const labelHashSet& sensitivityPatchIDs
        );

    // Destructor

       virtual ~adjointEikonalSolver() = default;


    // Member Functions

       //- Read dict if changed
       virtual bool readDict(const dictionary& dict);

       //- Accumulate source term
       void accumulateIntegrand(const scalar dt);

       //- Calculate the adjoint distance field
       void solve();

       //- Reset source term
       void reset();

       //- Return the sensitivity term depending on da
       boundaryVectorField& distanceSensitivities();

       //- Return the volume-based sensitivity term depending on da
       tmp<volTensorField> getFISensitivityTerm() const;

       //- Return sensitivity contribution to topology optimisation
       tmp<scalarField> topologySensitivities(const word& designVarsName) const;

       //-  Return the adjoint distance field
       const volScalarField& da();

       //- Return the distance field
       const volScalarField& d();

       //- Return the gradient of the eikonal equation
       tmp<volVectorField> gradEikonal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
