/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 PCOpt/NTUA
    Copyright (C) 2023 FOSS GP
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
    Foam::adjointNull

Description
    Dummy adjoint solver.
    Used to add the derivatives of geometric constraints which do not require
    the solution of adjoint PDEs

\*---------------------------------------------------------------------------*/

#ifndef adjointNull_H
#define adjointNull_H

#include "adjointSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class adjointNull Declaration
\*---------------------------------------------------------------------------*/

class adjointNull
:
    public adjointSolver
{

private:

    // Private Member Functions

        //- No copy construct
        adjointNull(const adjointNull&) = delete;

        //- No copy assignment
        void operator=(const adjointNull&) = delete;


protected:

    // Protected data

        //- Accumulate the sensitivities integrand before calculating them
        virtual void preCalculateSensitivities();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("null");


    // Constructors

        //- Construct from mesh and dictionary
        adjointNull
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointNull> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    //- Destructor
    virtual ~adjointNull() = default;


    // Member Functions

        // Access

            //- Return the type of simulation this solver pertains to
            virtual const word simulationType() const;

            //- Return the dimensions of the adjoint grid displacement variable
            virtual dimensionSet maDimensions() const;


        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Update primal based quantities related to the objective
            //- functions.
            //  Also writes the objective function values to files.
            virtual void updatePrimalBasedQuantities();

            // Functions related to the computation of sensitivity derivatives

                // Shape optimisation

                    //- Compute the multiplier for grad(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateGradDxDbMultiplier
                    (
                        volTensorField& gradDxDbMult,
                        const scalar dt
                    );

                    //- Compute the multiplier for div(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateDivDxDbMultiplier
                    (
                        autoPtr<scalarField>& divDxDbMult,
                        const scalar dt
                    );

                    //- Accumulate the multipliers of geometric quantities
                    //- defined at the boundary, usually through an objective
                    //- or constraint function
                    virtual void accumulateGeometryVariationsMultipliers
                    (
                        autoPtr<boundaryVectorField>& dSfdbMult,
                        autoPtr<boundaryVectorField>& dnfdbMult,
                        autoPtr<boundaryVectorField>& dxdbDirectMult,
                        autoPtr<pointBoundaryVectorField>& pointDxDirectDbMult,
                        const labelHashSet& sensitivityPatchIDs,
                        const scalar dt
                    );


                // Topology optimisation

                    //- Compute the multiplier of beta
                    virtual void topOSensMultiplier
                    (
                        scalarField& betaMult,
                        const word& designVariablesName,
                        const scalar dt
                    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
