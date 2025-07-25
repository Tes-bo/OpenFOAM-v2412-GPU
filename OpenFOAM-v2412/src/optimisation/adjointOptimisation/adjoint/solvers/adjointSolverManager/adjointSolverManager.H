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
    Foam::adjointSolverManager

Description
    Class for managing adjoint solvers, which may be more than one per
    operating point

SourceFiles
    adjointSolverManager.C

\*---------------------------------------------------------------------------*/

#ifndef adjointSolverManager_H
#define adjointSolverManager_H

#include "adjointSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class adjointSolverManager Declaration
\*---------------------------------------------------------------------------*/

class adjointSolverManager
:
    public regIOobject
{
private:

    // Private Member Functions

        //- No  copy construct
        adjointSolverManager(const adjointSolverManager&) = delete;

        //- No copy assignment
        void operator=(const adjointSolverManager&) = delete;


protected:

    // Protected Data

        fvMesh& mesh_;

        dictionary dict_;

        const word managerName_;

        const word managerType_;

        const word primalSolverName_;

        PtrList<adjointSolver> adjointSolvers_;

        labelList objectiveSolverIDs_;

        labelList oneSidedConstraintSolverIDs_;

        labelList doubleSidedConstraintSolverIDs_;

        scalar operatingPointWeight_;

        label nActiveAdjointSolvers_;

        autoPtr<designVariables>& designVars_;


public:

    TypeName("adjointSolverManager");


    // Constructors

        //- Construct from components
        adjointSolverManager
        (
            fvMesh& mesh,
            autoPtr<designVariables>& designVars,
            const word& managerType,
            const dictionary& dict,
            bool overrideUseSolverName
        );


    //- Destructor
    virtual ~adjointSolverManager() = default;


    // Member Functions

        bool readDict(const dictionary& dict);


        // Access

            //- Const access to the manager name
            const word& managerName() const;

            //- Const access to the primal solver name
            const word& primalSolverName() const;

            //- Const access to the construction dictionary
            const dictionary& dict() const;

            //- Const access to adjoint solvers
            const PtrList<adjointSolver>& adjointSolvers() const;

            //- Non-const access to adjoint solvers
            PtrList<adjointSolver>& adjointSolvers();

            //- Return the names of all adjointSolvers
            wordList adjointSolversNames() const;

            //- Const access to adjoint solvers
            scalar operatingPointWeight() const;

            //- Return number of active adjoint solvers, either corresponding
            //  to objectives or constraints
            label nActiveAdjointSolvers() const;

            //- Static function returning the number of active adjoint
            //  solvers reading dict
            static label nActiveAdjointSolvers(const dictionary& dict);

            //- Number of constraints
            //  Is the sum of one-sided and double-sided constraints, the
            //  latter multiplied by two
            label nConstraints() const;

            //- Number of adjoint solvers corresponding to one-sided
            //- constraints
            label nOneSidedConstraints() const;

            //- Number of adjoint solvers corresponding to double-sided
            //- constraints
            label nDoubleSidedConstraints() const;

            //- Number of adjoint solvers corresponding to objectives
            label nObjectives() const;

            //- Total number of adjoint solvers
            label nAdjointSolvers() const;


        // Evolution

            //- Update objective function-related values and solve adjoint
            //- equations
            void solveAdjointEquations();

            //- Aggregate sensitivities from various adjoint solvers
            tmp<scalarField> aggregateSensitivities();

            //- Get constraint sensitivities. One scalarField per constraint
            PtrList<scalarField> constraintSensitivities();

            //- Compute sensitivities for all adjoint solvers
            //- (both objective- and constraint-related ones)
            void computeAllSensitivities();

            //- Clear sensitivity fields from all adjoint solvers
            void clearSensitivities();

            //- Get objective value
            scalar objectiveValue();

            //- Get constraint values
            tmp<scalarField> constraintValues();

            //- Update fields related to primal solution.
            //  For instance, primal fields of adjoint turbulence models
            void updatePrimalBasedQuantities(const word& name);

            //- Whether the primal solver corresponding to the
            //- adjointSolverManager is the master one, in case of coupled
            //- solvers
            bool isMaster() const;


        // IO

            virtual bool writeData(Ostream&) const
            {
                return true;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
