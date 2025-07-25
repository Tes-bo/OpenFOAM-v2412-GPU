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
    Foam::objectives::objectiveFlowRatePartition

Description
    Split inlet flow rate to given percentages at the prescribed outlet patches.

SourceFiles
    objectiveFlowRatePartition.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveFlowRatePartition_H
#define objectiveFlowRatePartition_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                  Class objectiveFlowRatePartition Declaration
\*---------------------------------------------------------------------------*/

class objectiveFlowRatePartition
:
    public objectiveIncompressible
{
    // Private data

        labelList inletPatches_;
        labelList outletPatches_;
        scalarList targetFlowRateFraction_;
        scalarList currentFlowRateFraction_;
        scalar inletFlowRate_;
        scalarList flowRateDifference_;


public:

    //- Runtime type information
    TypeName("flowRatePartition");


    // Constructors

        //- from components
        objectiveFlowRatePartition
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveFlowRatePartition() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_boundarydJdv();

        //- Update values to be added to the adjoint outlet pressure
        virtual void update_boundarydJdvn();

        // Helper write functions

            //- Write any information that needs to go the header of the file
            //  (e.g. targets, directions, etc)
            virtual void addHeaderInfo() const;

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
