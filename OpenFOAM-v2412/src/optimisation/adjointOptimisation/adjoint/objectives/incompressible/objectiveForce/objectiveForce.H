/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2020, 2022 PCOpt/NTUA
    Copyright (C) 2013-2020, 2022 FOSS GP
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
    Foam::objectives::objectiveForce

Description

SourceFiles
    objectiveForce.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveForce_H
#define objectiveForce_H

#include "objectiveIncompressible.H"
#include "wallFvPatch.H"
#include "createZeroField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                        Class objectiveForce Declaration
\*---------------------------------------------------------------------------*/

class objectiveForce
:
    public objectiveIncompressible
{
protected:

    // Protected Data

        labelHashSet forcePatches_;
        vector forceDirection_;
        scalar Aref_;
        scalar rhoInf_;
        scalar UInf_;


public:

    //- Runtime type information
    TypeName("force");


    // Constructors

        //- Construct from components
        objectiveForce
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveForce() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint wall velocity
        virtual void update_boundarydJdp();

        //- Update delta(n dS)/delta b multiplier
        virtual void update_dSdbMultiplier();

        //- Update delta(x)/delta b multiplier
        virtual void update_dxdbMultiplier();

        //- Update dJ/dnut multiplier
        virtual void update_boundarydJdnut();

        //- Update dJ/dGradU multiplier
        virtual void update_boundarydJdGradU();

        //- Return denominator, without density
        virtual scalar denom() const;

        //- Return force direction
        virtual const vector& forceDirection() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
