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
    Foam::simple

Description
    Base class for solution control classes

\*---------------------------------------------------------------------------*/

#ifndef simple_H
#define simple_H

#include "incompressiblePrimalSolver.H"
#include "SIMPLEControl.H"
#include "IOMRFZoneList.H"
#include "objective.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class simple Declaration
\*---------------------------------------------------------------------------*/

class simple
:
    public incompressiblePrimalSolver
{
private:

    // Private Member Functions

        //- No copy construct
        simple(const simple&) = delete;

        //- No copy assignment
        void operator=(const simple&) = delete;


protected:

    // Protected data

        //- Solver control
        autoPtr<SIMPLEControl> solverControl_;

        //- Reference to incompressibleVars
        //  Used for convenience and to avoid repetitive dynamic_casts
        //  Same as getIncoVars()
        incompressibleVars& incoVars_;

        //- MRF zones
        IOMRFZoneList MRF_;

        //- Cumulative continuity error
        scalar cumulativeContErr_;

        //- List of objectives related to this primal solver
        UPtrList<objective> objectives_;

        //- Allow the execution of function objects in optimisation mode.
        //  Depending on the functions used and the setup, this might lead
        //  to a high amount of data written to the hard drive.
        //  WIP
        bool allowFunctionObjects_;


    //- Protected Member Functions

        //- Allocate incompressibleVars and return reference to be used for
        //- convenience in the rest of the class.
        incompressibleVars& allocateVars();

        //- In case variable names are different than the base ones,
        //- add extra schemes and relaxation factors to the appropriate dicts
        //  Note: Not supported for now
        void addExtraSchemes();

        //- Compute continuity errors
        void continuityErrors();


public:

    // Static Data Members

        //- Run-time type information
        TypeName("simple");


    // Constructors

        //- Construct from mesh and dictionary
        simple
        (
            fvMesh& mesh,
            const word& managerType,
            const dictionary& dict,
            const word& solverName
        );


    //- Destructor
    virtual ~simple() = default;


    // Member Functions

        virtual bool readDict(const dictionary& dict);

        // Evolution

            //- Execute one iteration of the solution algorithm
            virtual void solveIter();

            //- Steps to be executed before each main SIMPLE iteration
            virtual void preIter();

            //- The main SIMPLE iter
            virtual void mainIter();

            //- Steps to be executed before each main SIMPLE iteration
            virtual void postIter();

            //- Main control loop
            virtual void solve();

            //- Looper (advances iters, time step)
            virtual bool loop();

            //- Restore initial field values if necessary
            virtual void restoreInitValues();

            //- Functions to be called before loop
            virtual void preLoop();

            //- Functions to be called after loop
            virtual void postLoop();

            //- Write average iteration
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
