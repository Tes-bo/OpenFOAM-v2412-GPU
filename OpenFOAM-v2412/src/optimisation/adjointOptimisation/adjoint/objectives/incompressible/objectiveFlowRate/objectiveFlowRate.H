/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2022 PCOpt/NTUA
    Copyright (C) 2013-2022 FOSS GP
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
    Foam::objectives::objectiveFlowRate

Description
    Minimize/maximize flow rate through a given set of patches

SourceFiles
    objectiveFlowRate.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveFlowRate_H
#define objectiveFlowRate_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                      Class objectiveFlowRate Declaration
\*---------------------------------------------------------------------------*/

class objectiveFlowRate
:
    public objectiveIncompressible
{
    // Private data

        //- Patches used to define the objective function
        labelList patches_;

        //- Flow rate in each patch
        scalarList flowRates_;


public:

    //- Runtime type information
    TypeName("flowRate");


    // Constructors

        //- from components
        objectiveFlowRate
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveFlowRate() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_boundarydJdv();

        //- Update values to be added to the adjoint outlet pressure
        virtual void update_boundarydJdvn();

        // Helper write functions

            //- Write headers for additional columns
            virtual void addHeaderColumns() const;

            //- Write information to additional columns
            virtual void addColumnValues() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
