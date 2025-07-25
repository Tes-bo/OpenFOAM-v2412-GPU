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
    Foam::objectives::objectiveUniformityCellZone

Description
    Computes and minimizes (half) the variance of the velocity distribution in
    a given set of cellZones

SourceFiles
    objectiveUniformityCellZone.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveUniformityCellZone_H
#define objectiveUniformityCellZone_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                 Class objectiveUniformityCellZone Declaration
\*---------------------------------------------------------------------------*/

class objectiveUniformityCellZone
:
    public objectiveIncompressible
{
    // Private data

        //- Where to define the objective
        labelList zones_;

        //- Mean velocity value per objective zone
        vectorField UMean_;

        //- Velocity variance per objective zone
        scalarField UVar_;

        //- Volume integral per objective zone
        scalarField volZone_;


public:

    //- Runtime type information
    TypeName("uniformityCellZone");


    // Constructors

        //- from components
        objectiveUniformityCellZone
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveUniformityCellZone() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_dJdv();

        //- Update div( dx/db multiplier). Volume-based sensitivity term
        virtual void update_divDxDbMultiplier();

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
