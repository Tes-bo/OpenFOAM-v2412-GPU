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
    Foam::betaMaxValue

Description
    Read the betaMax value directly.

SourceFiles
    betaMaxValue.C

\*---------------------------------------------------------------------------*/

#ifndef betaMaxValue_H
#define betaMaxValue_H

#include "betaMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class betaMaxValue Declaration
\*---------------------------------------------------------------------------*/

class betaMaxValue
:
    public betaMax
{

private:

    // Private Member Functions

        //- No copy construct
        betaMaxValue(const betaMaxValue&) = delete;

        //- No copy assignment
        void operator=(const betaMaxValue&) = delete;


public:

    //- Runtime type information
    TypeName("value");


    // Constructors

        //- Construct from components
        betaMaxValue
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~betaMaxValue() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
