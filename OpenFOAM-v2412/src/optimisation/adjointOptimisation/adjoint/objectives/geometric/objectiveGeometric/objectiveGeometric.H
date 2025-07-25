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
    Foam::objectiveGeometric

Description
    Abstract base class for objective functions that contain only geometric
    quantities

SourceFiles
    objectiveGeometric.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveGeometric_H
#define objectiveGeometric_H

#include "objective.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class objectiveGeometric Declaration
\*---------------------------------------------------------------------------*/

class objectiveGeometric
:
    public objective
{

private:

    // Private Member Functions

        //- No copy construct
        objectiveGeometric(const objectiveGeometric&) = delete;

        //- No copy assignment
        void operator=(const objectiveGeometric&) = delete;


public:

    //- Runtime type information
    TypeName("geometric");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectiveGeometric,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict,
                const word& adjointSolverName,
                const word& primalSolverName
            ),
            (mesh, dict, adjointSolverName, primalSolverName)
        );


    // Constructors

        //- Construct from components
        objectiveGeometric
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objectiveGeometric> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveGeometric() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J() = 0;

        //- Update objective function derivatives
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
