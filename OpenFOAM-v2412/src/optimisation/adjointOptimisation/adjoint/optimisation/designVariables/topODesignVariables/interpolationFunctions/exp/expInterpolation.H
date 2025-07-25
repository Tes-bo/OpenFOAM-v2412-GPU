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
    Foam::expInterpolation

Description
    Computes the Brickmann penalization function for topology optimisation,
    relying on exp functions and proposed in

    Reference:
    \verbatim
        Wang, F., Lazarov, B., Sigmund, O.(2011).
        On projection methods, convergence and robust formulations in topology
        optimization.
        Structural and Multidisciplinary Optimization 43(6), 767–784.
        https://doi:10.1007/s00158-010-0602-y
    \endverbatim

SourceFiles
    expInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef expInterpolation_H
#define expInterpolation_H

#include "topOInterpolationFunction.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class expInterpolation Declaration
\*---------------------------------------------------------------------------*/

class expInterpolation
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
        expInterpolation(const expInterpolation&);

        //- No copy assignment
        void operator=(const expInterpolation&);


public:

    //- Runtime type information
    TypeName("exp");


    // Constructors

        //- Construct from components
        expInterpolation
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~expInterpolation() = default;


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
