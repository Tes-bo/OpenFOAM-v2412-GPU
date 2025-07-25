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
    Foam::objectiveManager

Description
    Class for managing objective functions.

SourceFiles
    objectiveManager.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveManager_H
#define objectiveManager_H

#include "fvMesh.H"
#include "objective.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class objectiveManager Declaration
\*---------------------------------------------------------------------------*/

class objectiveManager
:
    public regIOobject
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;
        const word adjointSolverName_;
        const word primalSolverName_;
        PtrList<objective> objectives_;
        autoPtr<OFstream> weightedObjectiveFile_;


private:

    // Private Member Functions

        //- No copy construct
        objectiveManager(const objectiveManager&) = delete;

        //- No copy assignment
        void operator=(const objectiveManager&) = delete;


public:

    TypeName("objectiveManager");

    // Constructors

        //- Construct from components
        objectiveManager
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveManager() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        //- Update objective function related values
        void updateNormalizationFactor();

        //- Update objective function related values
        void update();

        //- Update contributions to adjoint if true, otherwise return nulls
        void updateOrNullify();

        //- Increment integration times by the optimisation cycle time-span
        void incrementIntegrationTimes(const scalar timeSpan);

        //- Print to screen
        scalar print(bool negate = false);

        //- Should the objectives be written to file upon calling write()?
        void setWrite(const bool shouldWrite);
        //- Write objective function history
        virtual bool writeObjectives();

        //- Write objective function history
        virtual bool writeObjectives
        (
            const scalar weightedObjective,
            const bool valid = true
        );

        //- Call all functions required prior to the solution of the adjoint
        //- equations
        void updateAndWrite();

        //- Return reference to objective functions
        PtrList<objective>& getObjectiveFunctions();

        //- Return constant reference to objective functions
        const PtrList<objective>& getObjectiveFunctions() const;

        //- Return name of adjointSolverManager
        const word& adjointSolverManagerName() const;

        //- Return name of the adjointSolver
        const word& adjointSolverName() const;

        //- Return name of the primalSolver
        const word& primalSolverName() const;

        //- Check integration times for unsteady runs
        void checkIntegrationTimes() const;

        //- Add contribution to adjoint momentum PDEs
        virtual void addSource(fvVectorMatrix& matrix);

        //- Add contribution to a scalar adjoint PDEs
        virtual void addSource(fvScalarMatrix& matrix);


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
