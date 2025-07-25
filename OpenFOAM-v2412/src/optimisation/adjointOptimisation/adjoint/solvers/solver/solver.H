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
    Foam::solver

Description
    Base solver class

\*---------------------------------------------------------------------------*/

#ifndef solver_H
#define solver_H

#include "fvMesh.H"
#include "fvMatrix.H"
#include "localIOdictionary.H"
#include "variablesSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                            Class solver Declaration
\*---------------------------------------------------------------------------*/

class solver
:
    public localIOdictionary
{
protected:

    // Protected Data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- The optimisation type
        const word managerType_;

        //- Dictionary holding the solver info
        dictionary dict_;

        //- Solver name
        const word solverName_;

        //- Solve equations?
        bool active_;

        //- Is the solver the master one
        //  Used in cases where a solver needs to work in tandem with another
        bool isMaster_;

        //- Append the variables related to the solver with it name
        bool useSolverNameForFields_;

        //- Base variableSet pointer.
        //  To be allocated in derived classes
        autoPtr<variablesSet> vars_;


    // Protected Member Functions

        //- No copy construct
        solver(const solver&) = delete;

        //- No copy assignment
        void operator=(const solver&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solver");


    // Constructors

        //- Construct from mesh and dictionary
        solver
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& solverName
        );


    //- Destructor
    virtual ~solver() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);


        // Access

            //- Return the solver mesh
            inline const fvMesh& mesh() const;

            //- Return the solver name
            inline const word& solverName() const;

            //- Use solver name as a suffix to the involved fields
            inline bool useSolverNameForFields() const;

            //- Given a variable name, return a name that is possibly appended
            //- by the solverName, depending on useSolverNameForFields
            inline word extendedVariableName(const word& varName) const;

            //- Return state of solver
            inline bool active();

            //- Return the solver dictionary
            inline const dictionary& dict() const;

            //- Return constant reference to variableSet used by the solver
            inline const variablesSet& getVariablesSet() const;

            //- Return non-constant reference to variableSet used by the solver
            inline variablesSet& getVariablesSet();

            //- Return the manager type
            inline const word& managerType() const;

            //- Whether the solver is the master one, in case of coupled
            //- solvers
            inline bool isMaster() const;


        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter() = 0;

            //- Main control loop
            virtual void solve() = 0;

            //- Looper (advances iters, time step)
            virtual bool loop() = 0;

            //- Restore initial field values if necessary
            virtual void restoreInitValues();

            //- Functions to be called before loop
            virtual void preLoop();

            //- Functions to be called after loop
            virtual void postLoop();

            //- Main control loop.
            //  Gets a list of function pointers to be called at the end of
            //  each solver iteration
            template<class Type>
            void solveWithArgs
            (
                Type& type,
                List<void (Type::*)()>& funcs
            );

            //- Add topO fvOptions
            //  Does nothing in base
            virtual void addTopOFvOptions() const;


        // IO

            //- Required by regIOobject
            virtual bool writeData(Ostream&) const
            {
                return true;
            }

            //- Workaround for turbulent fields on multi-point runs
            virtual bool write(const bool valid = true) const
            {
                return false;
            }

            //- Workaround for turbulent fields on multi-point runs
            virtual bool writeNow() const
            {
                return false;
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "solverTemplates.C"
#endif
#include "solverI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
