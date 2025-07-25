/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2023 PCOpt/NTUA
    Copyright (C) 2013-2023 FOSS GP
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
    Foam::objectives::objectivePowerDissipation

Description
    Computes and minimizes the power dissipation within given cellZones.
    In the absence of significant viscous stresses on the "inlet" and "outlet"
    of the cellZones, this value is equal to the volume flow rate-weigthed
    total pressure losses (see also objectivePtLosses) within th cellZones

SourceFiles
    objectivePowerDissipation.C

\*---------------------------------------------------------------------------*/

#ifndef objectivePowerDissipation_H
#define objectivePowerDissipation_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                  Class objectivePowerDissipation Declaration
\*---------------------------------------------------------------------------*/

class objectivePowerDissipation
:
    public objectiveIncompressible
{
    // Private data

        //- Where to define the objective
        labelList zones_;


    // Private Member Functions

        //- Populate fieldNames
        void populateFieldNames();


public:

    //- Runtime type information
    TypeName("powerDissipation");


    // Constructors

        //- Construct from components
        objectivePowerDissipation
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectivePowerDissipation() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_dJdv();

        //- Update field to be added to the first adjoint turbulence model PDE
        virtual void update_dJdTMvar1();

        //- Update field to be added to the second adjoint turbulence model PDE
        virtual void update_dJdTMvar2();

        //- Update div(dx/db multiplier). Volume-based sensitivity term
        virtual void update_divDxDbMultiplier();

        //- Update grad(dx/db multiplier). Volume-based sensitivity term
        virtual void update_gradDxDbMultiplier();

        //- Contribution to field sensitivities
        virtual void update_dJdb();

        //- Add source terms to the adjoint turbulence model equations
        virtual void addSource(fvScalarMatrix& matrix);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
