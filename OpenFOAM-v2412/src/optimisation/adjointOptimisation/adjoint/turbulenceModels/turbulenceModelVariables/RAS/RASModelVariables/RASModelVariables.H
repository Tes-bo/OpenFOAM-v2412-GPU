/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019, 2022-2023 PCOpt/NTUA
    Copyright (C) 2013-2019, 2022-2023 FOSS GP
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
    Foam::incompressible::RASModelVariables

Description
    Abstract base class for objective functions. No point in making this
    runTime selectable since its children will have different constructors.

SourceFiles
    RASModelVariables.C

\*---------------------------------------------------------------------------*/

#ifndef RASModelVariables_H
#define RASModelVariables_H

#include "solverControl.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "runTimeSelectionTables.H"
#include "refPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{

/*---------------------------------------------------------------------------*\
                      Class RASModelVariables Declaration
\*---------------------------------------------------------------------------*/

class RASModelVariables
{
protected:

    // Protected Data

        const fvMesh& mesh_;
        const solverControl& solverControl_;

        // Base names of the turbulent fields
        word TMVar1BaseName_;
        word TMVar2BaseName_;
        word nutBaseName_;

        refPtr<volScalarField> TMVar1Ptr_;
        refPtr<volScalarField> TMVar2Ptr_;
        refPtr<volScalarField> nutPtr_;
        refPtr<volScalarField> distPtr_;

        // Initial values (demand-driven)
        // For finite differences and optimisation runs
        refPtr<volScalarField> TMVar1InitPtr_;
        refPtr<volScalarField> TMVar2InitPtr_;
        refPtr<volScalarField> nutInitPtr_;

        // Mean values (demand-driven)
        refPtr<volScalarField> TMVar1MeanPtr_;
        refPtr<volScalarField> TMVar2MeanPtr_;
        refPtr<volScalarField> nutMeanPtr_;


    // Protected Member functions

        virtual void allocateInitValues();
        virtual void allocateMeanFields();

        refPtr<volScalarField>
        cloneRefPtr(const refPtr<volScalarField>& obj) const;

        void copyAndRename(volScalarField& f1, volScalarField& f2);

        //- No copy assignment
        void operator=(const RASModelVariables&) = delete;


public:


    //- Runtime type information
    TypeName("RASModelVariables");

    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            RASModelVariables,
            dictionary,
            (
                const fvMesh& mesh,
                const solverControl& SolverControl
            ),
            (mesh, SolverControl)
        );


    // Constructors

        //- Construct from components
        RASModelVariables
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );

        //- Copy constructor
        //  Will allocate new fields (instead of referencing the ones in the
        //  turbulence model), so cannot be used directly to access the fields
        //  of the turbulence model.
        //  Mainly used for checkpointing in unsteady adjoint
        RASModelVariables
        (
            const RASModelVariables& rmv
        );

        //- Clone
        //  Will allocate new fields (instead of referencing the ones in the
        //  turbulence model), so cannot be used directly to access the fields
        //  of the turbulence model.
        //  Mainly used for checkpointing in unsteady adjoint
        autoPtr<RASModelVariables> clone() const;


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<RASModelVariables> New
        (
            const fvMesh& mesh,
            const solverControl& SolverControl
        );


    // Destructor

        // Destructor does nothing on base since depending on the case new
        // fields might be allocated
        // MUST be overloaded in inherited classes
        virtual ~RASModelVariables() = default;


    // Member Functions

        //- Turbulence field names
        inline const word& TMVar1BaseName() const;
        inline const word& TMVar2BaseName() const;
        inline const word& nutBaseName() const;

        //- Bools to identify which turbulent fields are present
        inline virtual bool hasTMVar1() const;
        inline virtual bool hasTMVar2() const;
        inline virtual bool hasNut() const;
        inline bool hasDist() const;

        //- Return references to turbulence fields
        //  will return the mean field if it exists, otherwise the
        //  instantaneous one
        inline const volScalarField& TMVar1() const;
        inline       volScalarField& TMVar1();
        inline const volScalarField& TMVar2() const;
        inline       volScalarField& TMVar2();
        inline const volScalarField& nutRef() const;
        inline       volScalarField& nutRef();
        inline tmp<volScalarField> nut() const;
        inline const volScalarField& d() const;
        inline       volScalarField& d();

        inline tmp<scalarField> TMVar1(const label patchi) const;
        inline tmp<scalarField> TMVar2(const label patchi) const;
        inline tmp<scalarField> nut(const label patchi) const;
        inline tmp<fvPatchScalarField> nutPatchField(const label patchi) const;

        //- Return references to instantaneous turbulence fields
        inline const volScalarField& TMVar1Inst() const;
        inline       volScalarField& TMVar1Inst();
        inline const volScalarField& TMVar2Inst() const;
        inline       volScalarField& TMVar2Inst();
        inline const volScalarField& nutRefInst() const;
        inline       volScalarField& nutRefInst();


        //- Return nut Jacobian wrt the TM vars
        virtual tmp<volScalarField> nutJacobianVar1
        (
            const singlePhaseTransportModel& laminarTransport
        ) const;

        virtual tmp<volScalarField> nutJacobianVar2
        (
            const singlePhaseTransportModel& laminarTransport
        ) const;

        //- Return the turbulence production term
        virtual tmp<volScalarField::Internal> G()
        {
            NotImplemented;
            return nullptr;
        }

        //- Restore turbulent fields to their initial values
        void restoreInitValues();

        //- Reset mean fields to zero
        void resetMeanFields();

        //- Compute mean fields on the fly
        virtual void computeMeanFields();

        //- Return stress tensor based on the mean flow variables
        tmp<volSymmTensorField> devReff
        (
            const singlePhaseTransportModel& laminarTransport,
            const volVectorField& U
        ) const;

        //- correct bounday conditions of turbulent fields
        virtual void correctBoundaryConditions
        (
            const incompressible::turbulenceModel& turbulence
        );

        //- Transfer turbulence fields from an another object
        //  Copies values since the ownership of the original fields is held
        //  by the turbulence model
        virtual void transfer(RASModelVariables& rmv);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RASModelVariablesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
