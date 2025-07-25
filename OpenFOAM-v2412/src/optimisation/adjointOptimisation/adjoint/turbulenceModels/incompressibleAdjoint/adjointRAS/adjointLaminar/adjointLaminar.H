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
    Foam::incompressibleAdjoint::adjointRASModels::adjointLaminar

Description
    Dummy turbulence model for a laminar incompressible flow. Can also
    be used when the "frozen turbulence" assumption is employed.

SourceFiles
    adjointLaminar.C

\*---------------------------------------------------------------------------*/

#ifndef adjointRasLaminar_H
#define adjointRasLaminar_H

#include "adjointRASModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressibleAdjoint
{
namespace adjointRASModels
{

/*---------------------------------------------------------------------------*\
                        Class adjointLaminar Declaration
\*---------------------------------------------------------------------------*/

class adjointLaminar
:
    public adjointRASModel
{
private:

    // Private Member Functions

        //- No copy construct
        adjointLaminar(const adjointLaminar&) = delete;

        //- No copy assignment
        void operator=(const adjointLaminar&) = delete;


public:

    //- Runtime type information
    TypeName("adjointLaminar");


    // Constructors

        //- Construct from components
        adjointLaminar
        (
            incompressibleVars& primalVars,
            incompressibleAdjointMeanFlowVars& adjointVars,
            objectiveManager& objManager,
            const word& adjointTurbulenceModelName
                = adjointTurbulenceModel::typeName,
            const word& modelName = typeName
        );


    //- Destructor
    virtual ~adjointLaminar() = default;


    // Member Functions

        //- Return the effective stress tensor, i.e. the adjointLaminar stress
        virtual tmp<volSymmTensorField> devReff() const;

        //- Return the effective stress tensor based on a given velocity field
        virtual tmp<volSymmTensorField> devReff
        (
            const volVectorField& U
        ) const;

        //- Return the diffusion term for the momentum equation
        virtual tmp<fvVectorMatrix> divDevReff(volVectorField& U) const;

        //- Source terms to the adjoint momentum equation due to the
        //- differentiation of the turbulence model
        virtual tmp<volVectorField> adjointMeanFlowSource();

        //- Returns zero field
        virtual const boundaryVectorField& adjointMomentumBCSource() const;

        //- Returns zero field
        virtual const boundaryVectorField& wallShapeSensitivities();

        //- Returns zero field
        virtual const boundaryVectorField& wallFloCoSensitivities();

        //- Returns zero field
        virtual tmp<volScalarField> distanceSensitivities();

        //- Returns zero field
        virtual tmp<volTensorField> FISensitivityTerm();

        //- Returns zero field
        virtual tmp<scalarField> topologySensitivities
        (
            const word& designVarsName
        ) const;

        //- Nullify all adjoint turbulence model fields and their old times
        virtual void nullify();

        //- Correct the primal viscosity field. Redundant?
        virtual void correct();

        //- Read adjointRASProperties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace adjointRASModels
} // End namespace incompressibleAdjoint
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
