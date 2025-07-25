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
    Foam::adjointSolver

Description
    Base class for adjoint solvers

\*---------------------------------------------------------------------------*/

#ifndef adjointSolver_H
#define adjointSolver_H

#include "fvMesh.H"
#include "Time.H"
#include "IOdictionary.H"
#include "solver.H"
#include "objectiveManager.H"
#include "primalSolver.H"
#include "adjointSensitivity.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class designVariables;

/*---------------------------------------------------------------------------*\
                        Class adjointSolver Declaration
\*---------------------------------------------------------------------------*/

class adjointSolver
:
    public solver
{
private:

    // Private Member Functions

        //- No copy construct
        adjointSolver(const adjointSolver&) = delete;

        //- No copy assignment
        void operator=(const adjointSolver&) = delete;


protected:

    // Protected data

        //- Name of primal solver
        const word primalSolverName_;

        //- Object to manage objective functions
        objectiveManager objectiveManager_;

        //- Sensitivities field
        tmp<scalarField> sensitivities_;

        //- Are sensitivities computed
        bool computeSensitivities_;

        //- Is the adjoint solver used to tackle a constraint
        bool isConstraint_;

        //- Is the adjoint solver used to tackle a double-sided constraint
        bool isDoubleSidedConstraint_;

        //- Sensitivity Derivatives engine
        autoPtr<adjointSensitivity> adjointSensitivity_;


    // Protected Member Functions

        //- Actions to be performed before calculating sensitivities
        //  Does noting in base
        virtual void preCalculateSensitivities()
        {}

        //- Allocate the sensitivity derivatives
        //  Since parts of the sensitivities depend on virtual functions
        //  implemented within derived classes, the actual allocation should
        //  happen there
        void allocateSensitivities();

        //- Return the dictionary corresponding to the design variables
        dictionary designVarsDict() const;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("adjointSolver");


    // Declare run-time constructor selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            adjointSolver,
            adjointSolver,
            (
                fvMesh& mesh,
                const word& managerType,
                const dictionary& dict,
                const word& primalSolverName,
                const word& solverName
            ),
            (mesh, managerType, dict, primalSolverName, solverName)
        );


    // Constructors

        //- Construct from mesh, dictionary, and primal solver name
        adjointSolver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointSolver> New
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& primalSolverName,
            const word& solverName
        );


    //- Destructor
    virtual ~adjointSolver() = default;


    // Member Functions

        // Access

            virtual bool readDict(const dictionary& dict);

            //- Return the primal solver name
            inline const word& primalSolverName() const;

            //- Return a const-reference to the primal solver
            //- corresponding to this adjoint solver
            inline const primalSolver& getPrimalSolver() const;

            //- Return a non const-reference to the primal solver
            //- corresponding to this adjoint solver
            inline primalSolver& getPrimalSolver();

            //- Return a const reference to the objective manager
            inline const objectiveManager& getObjectiveManager() const;

            //- Return a reference to the objective manager
            inline objectiveManager& getObjectiveManager();

            //- Is the solving referring to a constraint
            inline bool isConstraint();

            //- Is the solving referring to a double-sided constraint
            inline bool isDoubleSidedConstraint();

            //- Does the adjoint to an equation computing distances need to
            //- taken into consideration
            virtual bool includeDistance() const;

            //- Return the dimensions of the adjoint distance field
            virtual dimensionSet daDimensions() const;

            //- Return the dimensions of the adjoint grid displacement variable
            virtual dimensionSet maDimensions() const;

            //- Return the source the adjoint eikonal equation
            virtual tmp<volScalarField> adjointEikonalSource();

            //- Return the distance field, to be used in the solution of the
            //- adjoint eikonal PDE
            virtual tmp<volScalarField> yWall() const;


        // Evolution

            //- Compute sensitivities of the underlaying objectives
            virtual void computeObjectiveSensitivities
            (
                autoPtr<designVariables>& designVars
            );

            //- Grab a reference to the computed sensitivities
            virtual const scalarField& getObjectiveSensitivities
            (
                autoPtr<designVariables>& designVars
            );

            //- Clears the sensitivity field known by the adjoint solver
            virtual void clearSensitivities();

            //- Update primal based quantities, e.g. the primal fields
            //- in adjoint turbulence models
            //  Does nothing in the base
            virtual void updatePrimalBasedQuantities();

            //- Write the sensitivity derivatives
            virtual bool writeData(Ostream& os) const;

            // Functions related to the computation of sensitivity derivatives.
            // All functions get the field to accumulate their contribution on
            // as an argument and should be implemented by the derived classes

                // Shape optimisation

                    //- Compute the multiplier for grad(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateGradDxDbMultiplier
                    (
                        volTensorField& gradDxDbMult,
                        const scalar dt
                    )
                    {}

                    //- Compute the multiplier for div(dxdb)
                    //  Used in shape sensitivity derivatives, computed with
                    //  the FI and E-SI approaches
                    virtual void accumulateDivDxDbMultiplier
                    (
                        autoPtr<scalarField>& divDxDbMult,
                        const scalar dt
                    )
                    {}

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
                    )
                    {}

                    //- Contributions from boundary functions that inlcude
                    //- geometric aspects in them and change when the geometry
                    //- is displaced, e.g. rotationWallVelocity
                    virtual void accumulateBCSensitivityIntegrand
                    (
                        autoPtr<boundaryVectorField>& bcDxDbMult,
                        const labelHashSet& sensitivityPatchIDs,
                        const scalar dt
                    )
                    {}

                    //- Contributions from fvOptions that inlcude
                    //- geometric aspects in them and change when the geometry
                    //- is displaced, e.g. MRF
                    virtual void accumulateOptionsDxDbMultiplier
                    (
                        vectorField& optionsDxDbMult,
                        const scalar dt
                    )
                    {}


                // Topology optimisation

                    //- Compute the multiplier of beta
                    virtual void topOSensMultiplier
                    (
                        scalarField& betaMult,
                        const word& designVariablesName,
                        const scalar dt
                    )
                    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "adjointSolverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
