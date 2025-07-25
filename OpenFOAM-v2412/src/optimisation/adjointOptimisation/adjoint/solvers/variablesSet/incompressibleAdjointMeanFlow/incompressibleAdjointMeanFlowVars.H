/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2021 PCOpt/NTUA
    Copyright (C) 2013-2021 FOSS GP
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
    Foam::incompressibleAdjointMeanFlowVars

Description
    Manages the adjoint mean flow fields and their mean values

\*---------------------------------------------------------------------------*/

#ifndef incompressibleAdjointMeanFlowVars_H
#define incompressibleAdjointMeanFlowVars_H

#include "variablesSet.H"
#include "incompressibleVars.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class incompressibleAdjointMeanFlowVars Declaration
\*---------------------------------------------------------------------------*/

class incompressibleAdjointMeanFlowVars
:
    public variablesSet
{
private:

    // Private Member Functions

        //- No copy construct
        incompressibleAdjointMeanFlowVars
        (
            const incompressibleAdjointMeanFlowVars&
        ) = delete;

        //- No copy assignment
        void operator=(const incompressibleAdjointMeanFlowVars&) = delete;


protected:

    // Protected data

        //- Reference to the solverControl of the solver allocating the fields
        solverControl& solverControl_;

        //- Reference to primal variables
        incompressibleVars& primalVars_;

        //- Fields involved in the solution of the incompressible adjoint NS
        //- equations
        autoPtr<volScalarField> paPtr_;
        autoPtr<volVectorField> UaPtr_;
        autoPtr<surfaceScalarField> phiaPtr_;

        //- Mean Adjoint Fields. Actual averaging is done in the
        //- incompressibleAdjointVars class to take care of the mean adjoint
        //- turbulence variables
        autoPtr<volScalarField> paMeanPtr_;
        autoPtr<volVectorField> UaMeanPtr_;
        autoPtr<surfaceScalarField> phiaMeanPtr_;

        // Protected Member Functions

            //- Read fields and set turbulence
            void setFields();

            //- Read mean fields, if necessary
            void setMeanFields();


public:


    // Static Data Members

        //- Run-time type information
        TypeName("incompressibleAdjointMeanFlowVars");


    // Constructors

        //- Construct from mesh
        incompressibleAdjointMeanFlowVars
        (
            fvMesh& mesh,
            solverControl& SolverControl,
            incompressibleVars& primalVars
        );


    //- Destructor
    virtual ~incompressibleAdjointMeanFlowVars() = default;


    // Member Functions


        // Access

            // Const access to the primal vars related to the adjoint vars
            const incompressibleVars& primalVars() const;

            // Access to adjoint fields. Might be averaged or not depending on
            // the correspondign switch. Averaged fields are used to compute
            // an "average" sensitivity field

                //- Return const reference to pressure
                const volScalarField& pa() const;

                //- Return reference to pressure
                volScalarField& pa();

                //- Return const reference to velocity
                const volVectorField& Ua() const;

                //- Return reference to velocity
                volVectorField& Ua();

                //- Return const reference to volume flux
                const surfaceScalarField& phia() const;

                //- Return reference to volume flux
                surfaceScalarField& phia();


            // Access to instantaneous fields. Solvers and adjoint boundary
            // conditions should use these fields to execute a solver iteration

                //- Return const reference to pressure
                inline const volScalarField& paInst() const;

                //- Return reference to pressure
                inline volScalarField& paInst();

                //- Return const reference to velocity
                inline const volVectorField& UaInst() const;

                //- Return reference to velocity
                inline volVectorField& UaInst();

                //- Return const reference to volume flux
                inline const surfaceScalarField& phiaInst() const;

                //- Return reference to volume flux
                inline surfaceScalarField& phiaInst();

            //- Return computeMeanFields bool
            bool computeMeanFields() const;

            //- Return const reference to solverControl
            const solverControl& getSolverControl() const;

            //- Nullify all adjoint fields
            virtual void nullify();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "incompressibleAdjointMeanFlowVarsI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
