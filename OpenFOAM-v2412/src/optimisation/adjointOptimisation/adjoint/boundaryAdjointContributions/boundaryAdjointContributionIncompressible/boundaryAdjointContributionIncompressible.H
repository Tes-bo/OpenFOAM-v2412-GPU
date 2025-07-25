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
    Foam::boundaryAdjointContributionIncompressible

Description
    Contributions of objective function differentiation to adjoint
    boundary conditions for incompressible flows


SourceFiles
    boundaryAdjointContributionIncompressible.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryAdjointContributionIncompressible_H
#define boundaryAdjointContributionIncompressible_H

#include "boundaryAdjointContribution.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "objectiveManager.H"
#include "objectiveIncompressible.H"
#include "incompressibleAdjointSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
          Class boundaryAdjointContributionIncompressible Declaration
\*---------------------------------------------------------------------------*/

class boundaryAdjointContributionIncompressible
:
    public boundaryAdjointContribution
{
private:

    // Private Member Functions

        //- No copy construct
        boundaryAdjointContributionIncompressible
        (
            const boundaryAdjointContributionIncompressible&
        ) = delete;

        //- No copy assignment
        void operator=
        (
            const boundaryAdjointContributionIncompressible&
        ) = delete;


protected:

    // Protected data

        objectiveManager& objectiveManager_;

        const incompressibleVars& primalVars_;

        //- Note: getting a reference to the adjoint vars in the constructor of
        //- boundaryAdjointContributionIncompressible is dangerous since the
        //- autoPtr that holds them has not been completed yet. Instead, get
        //- a reference to the solver and grab the adjoint vars from there,
        //- when necessary
        const incompressibleAdjointSolver& adjointSolver_;


public:

    //- Runtime type information
    TypeName("incompressible");


    // Constructors

        //- Construct from components
        boundaryAdjointContributionIncompressible
        (
            const word& managerName,
            const word& adjointSolverName,
            const word& simulationType,
            const fvPatch& patch
        );


    //- Destructor
    virtual ~boundaryAdjointContributionIncompressible() = default;


    // Member Functions

//      tmp<vectorField> boundarydJtotdv();
        tmp<scalarField> pressureSource();
        tmp<vectorField> velocitySource();
        tmp<vectorField> tangentVelocitySource();
        tmp<vectorField> normalVelocitySource();
        tmp<scalarField> energySource();
        tmp<scalarField> adjointTMVariable1Source();
        tmp<scalarField> adjointTMVariable2Source();
        tmp<scalarField> dJdnut();
        tmp<tensorField> dJdGradU();

        tmp<scalarField> momentumDiffusion();
        tmp<scalarField> laminarDiffusivity();
        tmp<scalarField> thermalDiffusion();
        tmp<scalarField> wallDistance();

        tmp<scalarField> TMVariable1Diffusion();
        tmp<scalarField> TMVariable2Diffusion();
        tmp<scalarField> TMVariable1();
        tmp<scalarField> TMVariable2();

        const fvPatchVectorField& Ub() const;
        const fvPatchScalarField& pb() const;
        const fvsPatchScalarField& phib() const;
        tmp<fvPatchScalarField> turbulentDiffusivity() const;
        const fvPatchVectorField& Uab() const;
        const fvPatchScalarField& pab() const;
        const fvsPatchScalarField& phiab() const;

        const word primalSolverName() const;
        const word adjointSolverName() const;

        const incompressibleVars& primalVars()  const;
        const incompressibleAdjointVars& adjointVars() const;
        objectiveManager& getObjectiveManager();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
