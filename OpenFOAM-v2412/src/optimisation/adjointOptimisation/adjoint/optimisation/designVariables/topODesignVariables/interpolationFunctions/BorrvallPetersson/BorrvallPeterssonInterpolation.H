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
    Foam::BorrvallPeterssonInterpolation

Description
    Computes the Brickmann penalization function for topO optimisation,
    proposed by Borrvall and Petersson.

    Reference:
    \verbatim
        Borrvall, T., &  Petersson, J. (2002).
        TopO optimization of fluids in Stokes flow.
        International Journal for Numerical Methods in Fluids, 41(1), 77-107.
        https://doi.org/10.1002/fld.426
    \endverbatim

SourceFiles
    BorrvallPeterssonInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef BorrvallPeterssonInterpolation_H
#define BorrvallPeterssonInterpolation_H

#include "topOInterpolationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class BorrvallPeterssonInterpolation Declaration
\*---------------------------------------------------------------------------*/

class BorrvallPeterssonInterpolation
:
    public topOInterpolationFunction
{

private:

    // Private Data

        //- Steepness parameter
        //  High b values lead to more steep interpolations - higher contrasts.
        //  May vary in time
        autoPtr<Function1<scalar>> b_;


    // Private Member Functions

        //- No copy construct
        BorrvallPeterssonInterpolation
        (
            const BorrvallPeterssonInterpolation&
        ) = delete;

        //- No copy assignment
        void operator=(const BorrvallPeterssonInterpolation&) = delete;


public:

    //- Runtime type information
    TypeName("BorrvallPetersson");


    // Constructors

        //- Construct from components
        BorrvallPeterssonInterpolation
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~BorrvallPeterssonInterpolation() = default;


    // Member Functions

        //- Interpolate argument and write to result
        virtual void interpolate
        (
            const scalarField& arg,
            scalarField& res
        ) const;

        //- Return of function with respect to the argument field
        virtual tmp<scalarField> derivative(const scalarField& arg) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
