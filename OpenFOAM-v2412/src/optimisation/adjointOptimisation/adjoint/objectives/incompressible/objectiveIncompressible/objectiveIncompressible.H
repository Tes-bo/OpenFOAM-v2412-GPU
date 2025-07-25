/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::objectiveIncompressible

Description
    Abstract base class for objective functions in incompressible flows

SourceFiles
    objectiveIncompressible.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveIncompressible_H
#define objectiveIncompressible_H

#include "objective.H"
#include "incompressibleVars.H"
#include "adjointRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class objectiveIncompressible Declaration
\*---------------------------------------------------------------------------*/

class objectiveIncompressible
:
    public objective
{
protected:

    // Protected data

        const incompressibleVars& vars_;

        // Contribution to field adjoint equations
        // v,p,T and turbulence model variables
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<volVectorField> dJdvPtr_;
        autoPtr<volScalarField> dJdpPtr_;
        autoPtr<volScalarField> dJdTPtr_;

        //- First  turbulence model variable
        autoPtr<volScalarField> dJdTMvar1Ptr_;

        //- Second turbulence model variable
        autoPtr<volScalarField> dJdTMvar2Ptr_;

        // Contribution to adjoint boundary conditions
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<boundaryVectorField> bdJdvPtr_;

        //- Adjoint outlet pressure
        autoPtr<boundaryScalarField> bdJdvnPtr_;

        //- Adjoint outlet velocity
        autoPtr<boundaryVectorField> bdJdvtPtr_;

        //- Adjoint (intlet,wall) velocity
        autoPtr<boundaryVectorField> bdJdpPtr_;

        //- Adjoint outlet temperature
        autoPtr<boundaryScalarField> bdJdTPtr_;

        //- Adjoint outlet turbulence model var 1
        autoPtr<boundaryScalarField> bdJdTMvar1Ptr_;

        //- Adjoint outlet turbulence model var 2
        autoPtr<boundaryScalarField> bdJdTMvar2Ptr_;

        //- Jacobian wrt to nut
        autoPtr<boundaryScalarField> bdJdnutPtr_;

        //- Term multiplying gradU variations
        autoPtr<boundaryTensorField> bdJdGradUPtr_;


private:

    // Private Member Functions

        //- No copy construct
        objectiveIncompressible(const objectiveIncompressible&) = delete;

        //- No copy assignment
        void operator=(const objectiveIncompressible&) = delete;


public:

    //- Runtime type information
    TypeName("incompressible");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            objectiveIncompressible,
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
        objectiveIncompressible
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<objectiveIncompressible> New
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveIncompressible() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J() = 0;

        //- Normalize all fields allocated by the objective
        virtual void doNormalization();

        //- Contribution to field adjoint momentum eqs
        inline const volVectorField& dJdv();

        //- Contribution to field adjoint continuity eq
        inline const volScalarField& dJdp();

        //- Contribution to field adjoint energy eq
        inline const volScalarField& dJdT();

        //- Contribution to field adjoint turbulence model variable 1
        inline const volScalarField& dJdTMvar1();

        //- Contribution to field adjoint turbulence model variable 2
        inline const volScalarField& dJdTMvar2();

        //- Objective partial deriv wrt velocity for a specific patch
        const fvPatchVectorField& boundarydJdv(const label);

        //- Objective partial deriv wrt normal velocity for a specific patch
        inline const fvPatchScalarField& boundarydJdvn(const label);

        //- Objective partial deriv wrt tangent velocity for a specific patch
        inline const fvPatchVectorField& boundarydJdvt(const label);

        //- Objective partial deriv wrt pressure (times normal) for a specific
        //- patch
        inline const fvPatchVectorField& boundarydJdp(const label);

        //- Objective partial deriv wrt temperature for a specific patch
        inline const fvPatchScalarField& boundarydJdT(const label);

        //- Objective partial deriv wrt turbulence model var 1 for a specific
        //- patch
        inline const fvPatchScalarField& boundarydJdTMvar1(const label);

        //- Objective partial deriv wrt turbulence model var 2 for a specific
        //- patch
        inline const fvPatchScalarField& boundarydJdTMvar2(const label);

        //- Objective partial deriv wrt nut for a specific patch
        inline const fvPatchScalarField& boundarydJdnut(const label);

        //- Objective partial deriv wrt stress tensor
        inline const fvPatchTensorField& boundarydJdGradU(const label);

        //- Objective partial deriv wrt velocity for all patches
        inline const boundaryVectorField& boundarydJdv();

        //- Objective partial deriv wrt normal velocity for all patches
        inline const boundaryScalarField& boundarydJdvn();

        //- Objective partial deriv wrt tangent velocity for all patches
        inline const boundaryVectorField& boundarydJdvt();

        //- Objective partial deriv wrt pressure (times normal) for all patches
        inline const boundaryVectorField& boundarydJdp();

        //- Objective partial deriv wrt temperature for all patches
        inline const boundaryScalarField& boundarydJdT();

        //- Objective partial deriv wrt turbulence model var 1 for all patches
        inline const boundaryScalarField& boundarydJdTMvar1();

        //- Objective partial deriv wrt turbulence model var 2 for all patches
        inline const boundaryScalarField& boundarydJdTMvar2();

        //- Objective partial deriv wrt nut for all patches
        inline const boundaryScalarField& boundarydJdnut();

        //- Objective partial deriv wrt gradU
        inline const boundaryTensorField& boundarydJdGradU();

        //- Update objective function derivatives
        virtual void update();

        //- Update objective function derivatives
        virtual void nullify();

        // Funtions used in objectives including turbulent quantities

            //- Allocate fields related to the differentiation of turbulence
            //- models, if necessary
            void allocatedJdTurbulence();

            //- Check if cellZones provided include at least one cell
            void checkCellZonesSize(const labelList& zoneIDs) const;

            //- Compute dJdTMVar{1,2}
            void update_dJdTMvar
            (
                autoPtr<volScalarField>& dJdTMvarPtr,
                tmp<volScalarField>
                    (incompressibleAdjoint::adjointRASModel::*JacobianFunc)()
                    const,
                const volScalarField& JacobianMultiplier,
                const labelList& zones
            );


        //- Update vol and boundary fields and derivatives
        //  Do nothing in the base. The relevant ones should be overwritten
        //  in the child objective functions
        virtual void update_dJdv()
        {}

        virtual void update_dJdp()
        {}

        virtual void update_dJdT()
        {}

        virtual void update_dJdTMvar1()
        {}

        virtual void update_dJdTMvar2()
        {}

        virtual void update_dJdb()
        {}

        virtual void update_dJdbField()
        {}

        virtual void update_divDxDbMultiplier()
        {}

        virtual void update_gradDxDbMultiplier()
        {}

        virtual void update_boundarydJdv()
        {}

        virtual void update_boundarydJdvn()
        {}

        virtual void update_boundarydJdvt()
        {}

        virtual void update_boundarydJdp()
        {}

        virtual void update_boundarydJdT()
        {}

        virtual void update_boundarydJdTMvar1()
        {}

        virtual void update_boundarydJdTMvar2()
        {}

        virtual void update_boundarydJdnut()
        {}

        virtual void update_boundarydJdGradU()
        {}

        //- Vector sources can be given only to the adjoint momentum equations.
        //- Implemented in base objectiveIncompressible
        virtual void addSource(fvVectorMatrix& matrix);

        //- Scalar sources are more ambigious (adjoint pressure, turbulence
        //- model, energy, etc equations), so the equivalent functions should
        //- be overridden on an objective-basis
        virtual void addSource(fvScalarMatrix& matrix)
        {}


        //- Some objectives need to store some auxiliary values.
        //- If averaging is enabled, update these mean values here.
        //  By convention, the mean values (eg mean drag) refer to these flow
        //  values computed using the mean fields, rather than averaging the
        //  values themselves
        virtual void update_meanValues()
        {}

        //- Write objective function history
        virtual bool write(const bool valid = true) const;

        //- Inline functions for checking whether pointers are set or not
        inline bool hasdJdv() const noexcept;
        inline bool hasdJdp() const noexcept;
        inline bool hasdJdT() const noexcept;
        inline bool hasdJdTMVar1() const noexcept;
        inline bool hasdJdTMVar2() const noexcept;
        inline bool hasBoundarydJdv() const noexcept;
        inline bool hasBoundarydJdvn() const noexcept;
        inline bool hasBoundarydJdvt() const noexcept;
        inline bool hasBoundarydJdp() const noexcept;
        inline bool hasBoundarydJdT() const noexcept;
        inline bool hasBoundarydJdTMVar1() const noexcept;
        inline bool hasBoundarydJdTMVar2() const noexcept;
        inline bool hasBoundarydJdnut() const noexcept;
        inline bool hasBoundarydJdGradU() const noexcept;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "objectiveIncompressibleI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
