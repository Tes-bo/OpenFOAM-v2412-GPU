/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019, 2022 PCOpt/NTUA
    Copyright (C) 2013-2019, 2022 FOSS GP
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
    Foam::objectives::objectiveMoment

Description

SourceFiles
    objectiveMoment.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveMoment_H
#define objectiveMoment_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                       Class objectiveMoment Declaration
\*---------------------------------------------------------------------------*/

class objectiveMoment
:
    public objectiveIncompressible
{
    // Private data

        labelHashSet momentPatches_;
        vector momentDirection_;
        vector rotationCentre_;
        scalar Aref_;
        scalar lRef_;
        scalar rhoInf_;
        scalar UInf_;
        scalar invDenom_;

        // Store this in order to computed only once per objective call
        volSymmTensorField devReff_;


public:

    //- Runtime type information
    TypeName("moment");


    // Constructors

        //- from components
        objectiveMoment
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveMoment() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update mean drag and lift values
        virtual void update_meanValues();

        //- Update values to be added to the adjoint wall velocity
        virtual void update_boundarydJdp();

        //- Update delta(n dS)/delta b multiplier
        virtual void update_dSdbMultiplier();

        //- Update delta(x)/delta b multiplier
        virtual void update_dxdbMultiplier();

        //- Update delta(x)/delta b multiplier coming directly from the
        //- objective
        virtual void update_dxdbDirectMultiplier();

        //- Update dJ/dnut multiplier
        virtual void update_boundarydJdnut();

        //- Update dJ/dGradU multiplier
        /* WIP
        virtual void update_boundarydJdGradU();
        */
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectiveMoment
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
