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
    Foam::incompressible::adjointRASModels

Description
    Namespace for incompressible adjointRAS turbulence models.

Class
    Foam::incompressibleAdjoint::adjointRASModel

Description
    Abstract base class for incompressible turbulence models.

SourceFiles
    adjointRASModel.C

\*---------------------------------------------------------------------------*/

#ifndef adjointRASModel_H
#define adjointRASModel_H

#include "adjointTurbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "nearWallDist.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "incompressible/transportModel/transportModel.H"
#include "IOdictionary.H"
#include "Switch.H"
#include "bound.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "objectiveManager.H"
#include "boundaryFieldsFwd.H"
#include "createZeroField.H"
#include "solverControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{

/*---------------------------------------------------------------------------*\
                       Class adjointRASModel Declaration
\*---------------------------------------------------------------------------*/

class adjointRASModel
:
    public adjointTurbulenceModel,
    public IOdictionary
{
private:

    // Private Member Functions

        //- No copy construct
        adjointRASModel(const adjointRASModel&) = delete;

        //- No copy assignment
        void operator=(const adjointRASModel&) = delete;


protected:

    // Protected data

        //- Reference to the objectiveManager
        objectiveManager& objectiveManager_;

        //- Turbulence on/off flag
        Switch adjointTurbulence_;

        //- Flag to print the model coeffs at run-time
        Switch printCoeffs_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

        //- Near wall distance boundary field
        nearWallDist y_;

        //- Adjoint turbulence model variable 1
        autoPtr<volScalarField> adjointTMVariable1Ptr_;

        //- Adjoint turbulence model variable 2
        autoPtr<volScalarField> adjointTMVariable2Ptr_;

        //- Base names of the adjoint fields
        wordList adjointTMVariablesBaseNames_;

        //- Adjoint turbulence model variable 1, mean value
        autoPtr<volScalarField> adjointTMVariable1MeanPtr_;

        //- Adjoint turbulence model variable 2, mean value
        autoPtr<volScalarField> adjointTMVariable2MeanPtr_;

        //- Source to the adjoint momentum BC emerging
        //- from differentiating the turbulence model
        autoPtr<boundaryVectorField> adjMomentumBCSourcePtr_;

        //- Wall sensitivity term for shape optimisation
        autoPtr<boundaryVectorField> wallShapeSensitivitiesPtr_;

        //- Wall sensitivity term for flow control optimisation
        autoPtr<boundaryVectorField> wallFloCoSensitivitiesPtr_;

        //- Does the turbulence model include distances and should the
        //- adjoint to the distance field be computed
        bool includeDistance_;

        //- Has the primal solution changed?
        bool changedPrimalSolution_;


    // Protected Member Functions

        //- Print model coefficients
        virtual void printCoeffs();

        //- Set mean fields
        void setMeanFields();


public:

    //- Runtime type information
    TypeName("adjointRASModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            adjointRASModel,
            dictionary,
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
        adjointRASModel
        (
            const word& type,
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName =
                adjointTurbulenceModel::typeName
        );


    // Selectors

        //- Return a reference to the selected adjointRAS model
        static autoPtr<adjointRASModel> New
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName =
                adjointTurbulenceModel::typeName
        );


    //- Destructor
    virtual ~adjointRASModel() = default;


    // Member Functions

        //- Return the near wall distances
        const nearWallDist& y() const
        {
            return y_;
        }

        //- Const access to the coefficients dictionary
        const dictionary& coeffDict() const
        {
            return coeffDict_;
        }

        //- Const access to the primal solver name
        const word& primalSolverName() const
        {
            return primalVars_.solverName();
        }

        //- Const access to the adjoint solver name
        const word& adjointSolverName() const
        {
            return adjointVars_.solverName();
        }

        //- Return non-constant reference to adjoint turbulence model variable 1
        //  Will allocate and return a zero field in case it does not exist
        volScalarField& getAdjointTMVariable1Inst();

        //- Return non-constant reference to adjoint turbulence model variable 2
        //  Will allocate and return a zero field in case it does not exist
        volScalarField& getAdjointTMVariable2Inst();

        //- Return non-constant reference to adjoint turbulence model variable 1
        //  Will return the mean value if present,
        //  otherwise the instantaneous value
        volScalarField& getAdjointTMVariable1();

        //- Return non-constant reference to adjoint turbulence model variable 2
        //  Will return the mean value if present,
        //  otherwise the instantaneous value
        volScalarField& getAdjointTMVariable2();

        //- Return non-constant autoPtr to adjoint turbulence model variable 1
        autoPtr<volScalarField>& getAdjointTMVariable1InstPtr();

        //- Return non-constant autoPtr to adjoint turbulence model variable 2
        autoPtr<volScalarField>& getAdjointTMVariable2InstPtr();

        //- Return reference to the adjoint turbulence model variables base
        //- names
        const wordList& getAdjointTMVariablesBaseNames() const;

        //- Return the effective stress tensor including the laminar stress
        virtual tmp<volSymmTensorField> devReff() const = 0;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devReff
        (
            const volVectorField& U
        ) const = 0;

        //- Return the diffusion term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const = 0;

        //- Source terms to the adjoint momentum equation due to
        //- the differentiation of the turbulence model
        virtual tmp<volVectorField> adjointMeanFlowSource() = 0;

        //- Jacobian of nut wrt the first turbulence model variable
        //  Needed for objective functions that depend on nut. Defaults to zero
        virtual tmp<volScalarField> nutJacobianTMVar1() const;

        //- Jacobian of nut wrt the second turbulence model variable
        //  Needed for objective functions that depend on nut. Defaults to zero
        virtual tmp<volScalarField> nutJacobianTMVar2() const;

        //- Jacobian of nut wrt the flow velocity
        //  Assumes we want to get contributions of mult*d(nut)/dU
        //  Since the dependency of nut to U is usually through a differential
        //  operator, the term multiplying d(nut)/dU is passed as an argument
        //  to this function; the latter should then compute the
        //  contribution of the afforementioned term to the adjoint mean flow
        //  equations
        virtual tmp<volVectorField> nutJacobianU
        (
            tmp<volScalarField>& dNutdUMult
        ) const;

        //- Diffusion coefficient of the first primal and adjoint turbulence
        //- model equation. Needed for some adjoint BCs. Defaults to zero
        virtual tmp<scalarField> diffusionCoeffVar1(label patchI) const;

        //- Diffusion coefficient of the second primal and adjoint turbulence
        //- model equation. Needed for some adjoint BCs. Defaults to zero
        virtual tmp<scalarField> diffusionCoeffVar2(label patchI) const;

        //- Source for the outlet adjoint momentum BC coming from
        //- differentiating the turbulence model
        virtual const boundaryVectorField& adjointMomentumBCSource() const = 0;

        //- Sensitivity terms for shape optimisation, emerging from
        //- the turbulence model differentiation.
        //  Misses dxdb, to be added by the classes assembling the sensitivities
        virtual const boundaryVectorField& wallShapeSensitivities() = 0;

        //- Sensitivity terms for flow control, emerging from the
        //- turbulence model differentiation
        virtual const boundaryVectorField& wallFloCoSensitivities() = 0;

        //- Sensitivity terms resulting from the differentiation of the
        //- distance field. Misses dxdb, to be added by the classes
        //- assembling the sensitivities
        virtual tmp<volScalarField> distanceSensitivities() = 0;

        //- Term contributing to the computation of FI-based sensitivities
        //  Misses grad(dxdb), to be added by the assembling the sensitivities
        virtual tmp<volTensorField> FISensitivityTerm() = 0;

        //- Term contributing to the computation of topology optimisation
        //- sensitivities
        //  Misses betaMax, dBeta/dAlpha and the mesh volume, to be added
        //  during the assembly of the sensitivities
        virtual tmp<scalarField> topologySensitivities
        (
            const word& designVarsName
        ) const = 0;

        //- Solve the adjoint turbulence equations
        virtual void correct();

        //- Read adjointRASProperties dictionary
        virtual bool read();

        //- Set flag of changed primal solution to true
        void setChangedPrimalSolution();

        //- Restore field values to the initial ones
        void restoreInitValues();

        //- Reset mean fields to zero
        void resetMeanFields();

        //- Average adjoint fields on the fly
        void computeMeanFields();

        //- Should the adjoint to the eikonal equation be computed
        bool includeDistance() const;

        //- Nullify all adjoint turbulence model fields and their old times
        virtual void nullify() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif
