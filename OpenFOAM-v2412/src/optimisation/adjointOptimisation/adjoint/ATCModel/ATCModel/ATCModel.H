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
    Foam::ATCModel

Description
    Base class for selecting the adjoint transpose convection model.
    Inherits from regIOobject to add lookup functionality


SourceFiles
    ATCModel.C

\*---------------------------------------------------------------------------*/

#ifndef ATCModel_H
#define ATCModel_H

#include "regIOobject.H"
#include "autoPtr.H"
#include "zeroATCcells.H"
#include "incompressibleVars.H"
#include "incompressibleAdjointVars.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class ATCModel Declaration
\*---------------------------------------------------------------------------*/

class ATCModel
:
    public regIOobject
{
private:

    // Private Member Functions

        //- No copy construct
        ATCModel(const ATCModel&) = delete;

        //- No copy assignment
        void operator=(const ATCModel&) = delete;


protected:

    // Protected data

        const fvMesh& mesh_;
        const incompressibleVars& primalVars_;
        const incompressibleAdjointVars& adjointVars_;


        const dictionary& dict_;
        const scalar extraConvection_;
        const scalar extraDiffusion_;
        const label nSmooth_;
        const bool reconstructGradients_;
        word adjointSolverName_;
        autoPtr<zeroATCcells> zeroATCcells_;
        volScalarField ATClimiter_;
        volVectorField ATC_;


    // Protected functions

        //- Compute limiter based on the cells given by zeroATCcells
        void computeLimiter();

        //- Limit ATC field using ATClimiter_
        void smoothATC();


public:

    //- Runtime type information
    TypeName("ATCModel");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            ATCModel,
            dictionary,
            (
                const fvMesh& mesh,
                const incompressibleVars& primalVars,
                const incompressibleAdjointVars& adjointVars,
                const dictionary& dict
            ),
            (mesh, primalVars, adjointVars, dict)
        );


    // Constructors

        //- Construct from components
        ATCModel
        (
            const fvMesh& mesh,
            const incompressibleVars& primalVars,
            const incompressibleAdjointVars& adjointVars,
            const dictionary& dict
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<ATCModel> New
        (
            const fvMesh& mesh,
            const incompressibleVars& primalVars,
            const incompressibleAdjointVars& adjointVars,
            const dictionary& dict
        );


    //- Destructor
    virtual ~ATCModel() = default;


    // Member Functions

        //- Add ATC to the given matrix
        virtual void addATC(fvVectorMatrix& UaEqn) = 0;

        //- Get the list of cells on which to zero ATC
        const labelList& getZeroATCcells();

        //- Get the extra convection multiplier
        scalar getExtraConvectionMultiplier();

        //- Get the extra diffusion multiplier
        scalar getExtraDiffusionMultiplier();

        //- Get the list of cells on which to zero ATC
        const volScalarField& getLimiter() const;

        //- Compute limiter based on the prescribed cells and number of
        //- smoothing iterations
        static void computeLimiter
        (
            volScalarField& limiter,
            const labelList& smoothCells,
            const label nSmooth
        );

        //- Return a limiter based on given cells.
        //- For use with classes other than ATC which employ the same smoothing
        //- mechanism
        static tmp<volScalarField> createLimiter
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        //- Smooth an arbitrary field on a given list of cells
        template<class Type>
        void smoothFieldBasedOnCells
        (
            GeometricField<Type, fvPatchField, volMesh>& vf,
            const labelList& cells
        );

        //- Get the FI sensitivity derivatives term coming from the ATC
        virtual tmp<volTensorField> getFISensitivityTerm() const = 0;

        //- Dummy writeData function required from regIOobject
        virtual bool writeData(Ostream&) const;

        //- Update quantities related with the primal fields
        virtual void updatePrimalBasedQuantities();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "ATCModelTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
