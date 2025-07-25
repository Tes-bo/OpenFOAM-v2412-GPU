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
    Foam::solverControl

Description
    Base class for solver control classes

\*---------------------------------------------------------------------------*/

#ifndef solverControl_H
#define solverControl_H

#include "fvMesh.H"
#include "solver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class solverControl Declaration
\*---------------------------------------------------------------------------*/

class solverControl
{
protected:

    // Protected data

        //- Reference to the underlaying solver
        const solver& solver_;

        // Solution controls

            //- Whether to print the max magnitude during each solver iteration
            //  Useful for adjoint solvers
            bool printMaxMags_;


        // Evolution

            //- Current iteration index
            label iter_;

            //- Current averaging iteration index
            label averageIter_;

            //- Averaging start index
            label averageStartIter_;

            // Non run-time modifiable entries

                //- Whether to re-initialize the solution based on the initial
                //- time step. Makes sense only for optimisation or FD runs
                bool storeInitValues_;

                //- Do averaging
                bool average_;


    // Protected Member Functions

        //- Read controls from optimisationDict
        virtual bool read();


private:

        //- No copy construct
        solverControl(const solverControl&) = delete;

        //- No copy assignment
        void operator=(const solverControl&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solverControl");


    // Constructors

        //- Construct from solver
        solverControl(const solver& solver);


    //- Destructor
    virtual ~solverControl() = default;


    // Member Functions

        //- Read controls from optimisationDict

            //- Return the solver dictionary
            inline virtual const dictionary solverDict() const;

            //- Return the solutionControls dictionary
            inline virtual const dictionary solutionDict() const;

            //- Print max mags of solver fields
            inline bool printMaxMags() const;


        // Solution control

            //- Re-initialize
            inline bool storeInitValues() const;

            //- Return iteration index
            inline label iter() const;

            //- Return average iteration index reference
            inline label& averageIter();

            //- Return const average iteration index reference
            inline label averageIter() const;

            //- Return iteration index
            inline label averageStartIter() const;

            //- Whether or not to add fields of the current iteration to the
            //- average fields
            inline bool doAverageIter() const;

            //- Use averaged fields?
            //- For solving the adjoint equations or computing sensitivities
            //- based on averaged fields
            inline bool useAveragedFields() const;

            //- Whether averaging is enabled or not
            inline bool average() const;

        //- Return reference to the underlaying solver
        inline const solver& getSolver() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#include "solverControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
