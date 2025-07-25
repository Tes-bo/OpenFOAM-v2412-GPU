/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 PCOpt/NTUA
    Copyright (C) 2020-2023 FOSS GP
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
    Foam::betaMaxStepRamp

Description
    Computes the betaMax using a stepRamp function which updates itself
    based on the current time

\*---------------------------------------------------------------------------*/

#ifndef betaMaxStepRamp_H
#define betaMaxStepRamp_H

#include "betaMax.H"
#include "stepRamp.H"
#include "Scale.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class betaMaxStepRamp Declaration
\*---------------------------------------------------------------------------*/

class betaMaxStepRamp
:
    public betaMax
{
private:

    // Private Member Functions

        //- No copy construct
        betaMaxStepRamp(const betaMaxStepRamp&) = delete;

        //- No copy assignment
        void operator=(const betaMaxStepRamp&) = delete;


protected:

    // Protected Data Members

        //- Minimum betaMax value (to be used at the first opt cycles)
        scalar betaMin_;

        //- The step ramp function
        autoPtr<Function1<scalar>> funcPtr_;


public:

    //- Runtime type information
    TypeName("stepRamp");


    // Constructors

        //- Construct from components
        betaMaxStepRamp
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~betaMaxStepRamp() = default;


    // Member Functions

        //- Get value
        virtual scalar value() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
