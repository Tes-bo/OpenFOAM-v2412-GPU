/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::objectives::objectiveTopOSolidVolume

Description
    Objective quantifying the difference between the volume occupied by solid
    in topology optimisation and a target percentage; the latter can change
    throughout the optimisation cycles through a Function1.

SourceFiles
    objectiveTopOSolidVolume.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveTopOSolidVolume_H
#define objectiveTopOSolidVolume_H

#include "objectiveGeometric.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                   Class objectiveTopOSolidVolume Declaration
\*---------------------------------------------------------------------------*/

class objectiveTopOSolidVolume
:
    public objectiveGeometric
{
    // Private data

        autoPtr<Function1<scalar>> targetPercentage_;

        bool percentInDenom_;


public:

    //- Runtime type information
    TypeName("topOSolidVolume");


    // Constructors

        //- from components
        objectiveTopOSolidVolume
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveTopOSolidVolume() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Contribution to field sensitivities
        virtual void update_dJdb();

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
