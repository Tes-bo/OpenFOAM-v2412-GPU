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


Namespace
    Foam::incompressibleAdjoint

Description
    Namespace for incompressible adjoint turbulence models.

Class
    Foam::incompressibleAdjoint::adjointTurbulenceModel

Description
    Abstract base class for incompressible adjoint turbulence models
    (RAS, LES and laminar).

SourceFiles
    adjointTurbulenceModel.C
    newTurbulenceModel.C

\*---------------------------------------------------------------------------*/

#ifndef adjointTurbulenceModel_H
#define adjointTurbulenceModel_H

#include "incompressibleVars.H"
#include "incompressibleAdjointMeanFlowVars.H"
#include "objectiveManager.H"
#include "Time.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class fvMesh;

namespace incompressibleAdjoint
{

/*---------------------------------------------------------------------------*\
                    Class adjointTurbulenceModel Declaration
\*---------------------------------------------------------------------------*/

class adjointTurbulenceModel
:
    public regIOobject
{
private:

    // Private Member Functions

        //- No copy construct
        adjointTurbulenceModel(const adjointTurbulenceModel&) = delete;

        //- No copy assignment
        void operator=(const adjointTurbulenceModel&) = delete;


protected:

    // Protected data

        incompressibleVars& primalVars_;
        incompressibleAdjointMeanFlowVars& adjointVars_;
        const Time& runTime_;
        const fvMesh& mesh_;


public:

    //- Runtime type information
    TypeName("adjointTurbulenceModel");


    // Declare run-time New selection table

        declareRunTimeNewSelectionTable
        (
            autoPtr,
            adjointTurbulenceModel,
            adjointTurbulenceModel,
            (
                incompressibleVars& primalVars,
                incompressibleAdjointMeanFlowVars& adjointVars,
                objectiveManager& objManager,
                const word& adjointTurbulenceModelName
            ),
            (
                primalVars,
                adjointVars,
                objManager,
                adjointTurbulenceModelName
            )
        );


    // Constructors

        //- Construct from components
        adjointTurbulenceModel
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName = typeName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<adjointTurbulenceModel> New
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName = typeName
        );


    //- Destructor
    virtual ~adjointTurbulenceModel() = default;


    // Member Functions

        //- Return the laminar viscosity
        inline tmp<volScalarField> nu() const
        {
            return primalVars_.laminarTransport().nu();
        }

        //- Return the turbulence viscosity
        virtual const tmp<volScalarField> nut() const
        {
            return primalVars_.RASModelVariables()().nut();
        }

        //- Return the effective viscosity
        virtual tmp<volScalarField> nuEff() const
        {
            // Go through RASModelVariables::nutRef in order to obtain
            // the mean field, if present
            const singlePhaseTransportModel& lamTrans =
                primalVars_.laminarTransport();
            const autoPtr<incompressible::RASModelVariables>&
                turbVars = primalVars_.RASModelVariables();

            return volScalarField::New
            (
                "nuEff",
                IOobject::NO_REGISTER,
                lamTrans.nu() + turbVars().nut()
            );
        }

        //- Return the effective viscosity on a given patch
        virtual tmp<scalarField> nuEff(const label patchI) const
        {
            // Go through RASModelVariables::nutRef in order to obtain
            // the mean field, if present
            const singlePhaseTransportModel& lamTrans =
                primalVars_.laminarTransport();
            const autoPtr<incompressible::RASModelVariables>&
                turbVars = primalVars_.RASModelVariables();

            return (lamTrans.nu(patchI) + turbVars().nut(patchI));
        }

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devReff
        (
            const volVectorField& U
        ) const = 0;

        //- Return the diffusion term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

        //- Source term added to the adjoint mean flow due to the
        //  differentiation of the turbulence model
        virtual tmp<volVectorField> adjointMeanFlowSource() = 0;

        //- Solve the adjoint turbulence equations
        virtual void correct() = 0;

        //- Read adjointLESProperties or adjointRASProperties dictionary
        virtual bool read() = 0;

        //- Default dummy write function
        virtual bool writeData(Ostream&) const
        {
            return true;
        }

        //- Nullify all adjoint turbulence model fields and their old times
        virtual void nullify() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressibleAdjoint
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
