/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Interpolation CFD Toolbox
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
    Foam::SIMPInterpolation

Description
    Computes the Brickmann penalization function for topology optimisation,
    using the SIMP approach proposed by Bendsøe

    Reference:
    \verbatim
        Bendsøe, M. P. (1989)
        Optimal shape design as a material distribution problem.
        Structural optimization, 1(4), 193–202.
        https://doi.org/10.1007/BF01650949
    \endverbatim

InterpolationFiles
    SIMPInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef SIMPInterpolation_H
#define SIMPInterpolation_H

#include "topOInterpolationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class SIMPInterpolation Declaration
\*---------------------------------------------------------------------------*/

class SIMPInterpolation
:
    public topOInterpolationFunction
{
    // Private Data

        //- Steepness parameter
        //  High b values lead to more steep interpolations - higher contrasts.
        //  May vary in time
        autoPtr<Function1<scalar>> b_;

private:

    // Private Member Functions

        //- No copy construct
        SIMPInterpolation(const SIMPInterpolation&);

        //- No copy assignment
        void operator=(const SIMPInterpolation&);


public:

    //- Runtime type information
    TypeName("SIMP");


    // Constructors

        //- Construct from components
        SIMPInterpolation
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~SIMPInterpolation() = default;


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
