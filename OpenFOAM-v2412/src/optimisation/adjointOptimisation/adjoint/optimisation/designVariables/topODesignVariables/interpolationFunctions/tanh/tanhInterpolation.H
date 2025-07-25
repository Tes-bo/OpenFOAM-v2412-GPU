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
    Foam::tanhInterpolation

Description
    Computes the Brickmann penalization function for topology optimisation,
    relying on tanh functions and proposed in

    Reference:
    \verbatim
        Wang, F., Lazarov, B., Sigmund, O.(2011).
        On projection methods, convergence and robust formulations in topO
        optimization.
        Structural and Multidisciplinary Optimization 43(6), 767–784.
        https://doi:10.1007/s00158-010-0602-y
    \endverbatim

SourceFiles
    tanhInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef tanhInterpolation_H
#define tanhInterpolation_H

#include "topOInterpolationFunction.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class tanhInterpolation Declaration
\*---------------------------------------------------------------------------*/

class tanhInterpolation
:
    public topOInterpolationFunction
{

private:

    // Private Data

        //- Steepness parameter
        //  High b values lead to more steep interpolations - higher contrasts.
        //  May vary in time
        autoPtr<Function1<scalar>> b_;

        //- Length scale parameter
        //  Default to 0.5
        scalar eta_;


    // Private Member Functions

        //- No copy construct
        tanhInterpolation(const tanhInterpolation&);

        //- No copy assignment
        void operator=(const tanhInterpolation&);


public:

    //- Runtime type information
    TypeName("tanh");


    // Constructors

        //- Construct from components
        tanhInterpolation
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~tanhInterpolation() = default;


    // Member Functions

        //- Interpolate argument and write to result
        virtual void interpolate
        (
            const scalarField& arg,
            scalarField& res
        ) const;

        //- Return of function with respect to the argument field
        virtual tmp<scalarField> derivative(const scalarField& arg) const;

        //- Set the parameter determining length scale
        virtual void setLengthScaleParam(const scalar lengthScale);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
