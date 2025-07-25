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
    Foam::objectives::objectiveNutSqr

Description
    Objective qualitatively quantifying noise through the integral of the
    squared turbulent viscosity over specified cellZones. Requires the adjoint
    to the turbulence model to be incorporated into the optimisation loop.

    Reference:
    \verbatim
        Objective function initially presented in

            Papoutsis-Kiachagias, E. M., Magoulas, N., Mueller, J., Othmer, C.,
            & Giannakoglou, K. C. (2015).
            Noise reduction in car aerodynamics using a surrogate objective
            function and the continuous adjoint method with wall functions.
            Computers & Fluids, 122(20), 223-232.
            https://doi.org/10.1016/j.compfluid.2015.09.002
    \endverbatim

SourceFiles
    objectiveNutSqrNoise.C

\*---------------------------------------------------------------------------*/

#ifndef objectiveNutSqr_H
#define objectiveNutSqr_H

#include "objectiveIncompressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace objectives
{

/*---------------------------------------------------------------------------*\
                       Class objectiveNutSqr Declaration
\*---------------------------------------------------------------------------*/

class objectiveNutSqr
:
    public objectiveIncompressible
{
    // Private data

        //- Where to define the objective
        labelList zones_;

        //- List with the names of the adjoint turbulence model fields
        //  This is kept separately from fieldNames since the latter may
        //  or may not include Ua, depending on whether nut is a function
        //  of U.  This makes deciding on whether to add or not sources
        //  to a given fvScalarMatrix tricky, hence the utilisation of
        //  adjointTurbulenceNames_
        wordList adjointTurbulenceNames_;


    // Private Member Functions

        //- Populate fieldNames
        void populateFieldNames();


public:

    //- Runtime type information
    TypeName("nutSqr");


    // Constructors

        //- Construct from components
        objectiveNutSqr
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName,
            const word& primalSolverName
        );


    //- Destructor
    virtual ~objectiveNutSqr() = default;


    // Member Functions

        //- Return the objective function value
        virtual scalar J();

        //- Update values to be added to the adjoint outlet velocity
        virtual void update_dJdv();

        //- Update field to be added to the first adjoint turbulence model PDE
        virtual void update_dJdTMvar1();

        //- Update field to be added to the second adjoint turbulence model PDE
        virtual void update_dJdTMvar2();

        //- Update field to be added to be added to volume-based
        //- sensitivity derivatives, emerging from delta ( dV ) / delta b
        virtual void update_divDxDbMultiplier();

        //- Add source terms to the adjoint turbulence model equations
        virtual void addSource(fvScalarMatrix& matrix);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace objectives
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
