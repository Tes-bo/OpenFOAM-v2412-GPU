/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020 PCOpt/NTUA
    Copyright (C) 2013-2020 FOSS GP
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::objectives::objectivePtLosses

Description

SourceFiles
    objectivePtLosses.C

\*---------------------------------------------------------------------------*/

#ifndef objectivePtLosses_H
#define objectivePtLosses_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                      Class objectivePtLosses Declaration
\*---------------------------------------------------------------------------*/

class objectivePtLosses
:
    public objectiveIncompressible
{
    // Private data

        //- Patches on which to compute losses
        labelList patches_;

        //- Total pressure per patch on patches_
        scalarField patchPt_;


    // Private Member Functions

        //- Return the objectiveReportPatches
        void initialize();


public:

    //- Runtime type information
    TypeName("PtLosses");


    // Constructors

        //- from components
        objectivePtLosses
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectivePtLosses() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint inlet velocity
        virtual void update_boundarydJdp();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_boundarydJdv();

        //- Update values to be added to the adjoint outlet pressure
        virtual void update_boundarydJdvn();

        //- Update values to be added to the adjoint outlet tangential velocity
        virtual void update_boundarydJdvt();

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
