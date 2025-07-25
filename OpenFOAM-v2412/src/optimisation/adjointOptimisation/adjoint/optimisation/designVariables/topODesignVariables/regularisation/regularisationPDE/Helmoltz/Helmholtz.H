/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2023 PCOpt/NTUA
    Copyright (C) 2021-2023 FOSS GP
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
    Foam::Helmholtz

Description
    A regulatisation PDE based on a Helmholtz filter


SourceFiles
    Helmholtz.C

\*---------------------------------------------------------------------------*/

#ifndef Helmholtz_H
#define Helmholtz_H

#include "regularisationPDE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class Helmholtz Declaration
\*---------------------------------------------------------------------------*/

class Helmholtz
:
    public regularisationPDE
{

private:

    // Private Member Functions

        //- No copy construct
        Helmholtz(const Helmholtz&) = delete;

        //- No copy assignment
        void operator=(const Helmholtz&) = delete;


protected:

    // Protected Data

        //- Solve the regularisationPDE only on a subset mesh made of the
        //- active cell zones
        bool solveOnActiveCells_;

        //- Fixed value at wall boundaries. Defaults to 1
        scalar wallValue_;


    // Protected Member Functions

        //- Solve the regulatisation equation
        //  The mesh used is the one of aTilda
        void solveEqn
        (
            const volScalarField& aTilda,
            const scalarField& source,
            scalarField& result,
            const bool isTopoField,
            const regularisationRadius& radius,
            const scalar minSetValue = Zero,
            const bool fixATildaValues = true
        );


public:

    //- Runtime type information
    TypeName("Helmholtz");


    // Constructors

        //- Construct from components
        Helmholtz
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const topOZones& zones
        );


    //- Destructor
    virtual ~Helmholtz() = default;


    // Member Functions

        //- Regularise field (or sensitivities) using a Laplace PDE
        virtual void regularise
        (
            const volScalarField& aTilda,
            const scalarField& source,
            scalarField& result,
            const bool isTopoField,
            const regularisationRadius& radius,
            const scalar minSetValue = Zero,
            const bool fixATildaValues = true
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
