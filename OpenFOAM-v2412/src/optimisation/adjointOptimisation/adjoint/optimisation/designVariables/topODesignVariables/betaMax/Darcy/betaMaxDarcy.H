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
    Foam::betaMaxDarcy

Description
    Compute betaMax through the definition of the Darcy number, quantifying
    the viscous-to-porous forces ratio
        Da = nu/betaMax/L/L
    where nu is the bulk viscosity and L is a characteristic length.
    The latter is either supplied directly or computed as the
    - 2D: area of the inlet patches divided by the span in the empty direction
    - 3D: the hydraulic diameter of the inlet patches

    Reference:
    \verbatim
        Olesen, L. H., Okkels, F., & Bruus H. (2006).
        A high-level programming-language implementation of topology
        optimization applied to steady-state Navier–Stokes flow.
        International Journal for Numerical Methods in Engineering,
        65(7), 975-1001.
        https://doi.org/10.1002/nme.1468
    \endverbatim

SourceFiles
    betaMaxDarcy.C

\*---------------------------------------------------------------------------*/

#ifndef betaMaxDarcy_H
#define betaMaxDarcy_H

#include "betaMax.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class betaMaxDarcy Declaration
\*---------------------------------------------------------------------------*/

class betaMaxDarcy
:
    public betaMax
{

private:

    // Private Member Functions

        //- No copy construct
        betaMaxDarcy(const betaMaxDarcy&) = delete;

        //- No copy assignment
        void operator=(const betaMaxDarcy&) = delete;


protected:

    // Protected Data

        //- The Darcy number expressing the ratio of viscous to porous forces
        scalar DarcyNumber_;

        //- Characteristic length of the case
        //  Either supplied directly or computed as the hydraulic diameter of
        //  the inlet
        scalar length_;


public:

    //- Runtime type information
    TypeName("Darcy");


    // Constructors

        //- Construct from components
        betaMaxDarcy
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~betaMaxDarcy() = default;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
