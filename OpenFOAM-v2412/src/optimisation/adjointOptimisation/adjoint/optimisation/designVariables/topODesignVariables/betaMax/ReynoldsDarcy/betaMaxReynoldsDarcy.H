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
    Foam::betaMaxReynoldsDarcy

Description
    Compute betaMax through the products of the Reynolds and Darcy numbers,
    quantifying the momentum-to-porous forces ratio
        ReDa = Uref/betaMax/L
    where Uref is a reference velocity and L is a characteristic length.
    The latter is either supplied directly or computed as the
    - 2D: area of the inlet patches divided by the span in the empty direction
    - 3D: the hydraulic diameter of the inlet patches

SourceFiles
    betaMaxReynoldsDarcy.C

\*---------------------------------------------------------------------------*/

#ifndef betaMaxReynoldsDarcy_H
#define betaMaxReynoldsDarcy_H

#include "betaMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class betaMaxReynoldsDarcy Declaration
\*---------------------------------------------------------------------------*/

class betaMaxReynoldsDarcy
:
    public betaMax
{

private:

    // Private Member Functions

        //- No copy construct
        betaMaxReynoldsDarcy(const betaMaxReynoldsDarcy&) = delete;

        //- No copy assignment
        void operator=(const betaMaxReynoldsDarcy&) = delete;


protected:

    // Protected Data

        //- The Darcy number expressing the ratio of viscous to porous forces
        scalar ReynoldsDarcyNumber_;

        //- Characteristic length of the case
        //  Either supplied directly or computed as the hydraulic diameter of
        //  the inlet
        scalar length_;

        //- Reference velocity
        scalar Uref_;


public:

    //- Runtime type information
    TypeName("ReynoldsDarcy");


    // Constructors

        //- Construct from components
        betaMaxReynoldsDarcy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~betaMaxReynoldsDarcy() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
