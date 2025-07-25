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
    Foam::incompressibleAdjointVars

Description
    Class including all adjoint fields for incompressible flows

\*---------------------------------------------------------------------------*/

#ifndef incompressibleAdjointVars_H
#define incompressibleAdjointVars_H

#include "incompressibleAdjointMeanFlowVars.H"
#include "objectiveManager.H"
#include "adjointRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class incompressibleAdjointVars Declaration
\*---------------------------------------------------------------------------*/

class incompressibleAdjointVars
:
    public incompressibleAdjointMeanFlowVars
{
protected:

    // Protected data

        //- Reference to the objectiveManager
        objectiveManager& objectiveManager_;

        //- Adjoint to the turbulence model
        autoPtr<incompressibleAdjoint::adjointRASModel> adjointTurbulence_;


    // Protected Member Functions

        //- No copy construct
        incompressibleAdjointVars(const incompressibleAdjointVars&) = delete;

        //- No copy assignment
        void operator=(const incompressibleAdjointVars&) = delete;


public:


    // Static Data Members

        //- Run-time type information
        TypeName("incompressibleAdjointVars");


    // Constructors

        //- Construct from mesh
        incompressibleAdjointVars
        (
            fvMesh& mesh,
            solverControl& SolverControl,
            objectiveManager& objManager,
            incompressibleVars& primalVars
        );


    //- Destructor
    virtual ~incompressibleAdjointVars() = default;


    // Member Functions

        // Access

            //- Return const reference to the adjointRASModel
            inline const autoPtr<incompressibleAdjoint::adjointRASModel>&
                adjointTurbulence() const;

            //- Return non-const reference to the adjointRASModel
            inline autoPtr<incompressibleAdjoint::adjointRASModel>&
                adjointTurbulence();

            //- Restore field values to the initial ones
            void restoreInitValues();

            //- Reset mean fields to zero
            void resetMeanFields();

            //- Compute mean fields on the fly
            void computeMeanFields();

            //- Nullify all adjoint fields
            virtual void nullify();

            //- Update primal based quantities of the adjoint boundary
            //  conditions
            virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "incompressibleAdjointVarsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
