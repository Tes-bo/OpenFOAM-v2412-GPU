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
    Foam::linearInterpolation

Description
    Implements a linear (with respect to beta) term for topology optimisation

InterpolationFiles
    linearInterpolation.C

\*---------------------------------------------------------------------------*/

#ifndef linearInterpolation_H
#define linearInterpolation_H

#include "topOInterpolationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class linearInterpolation Declaration
\*---------------------------------------------------------------------------*/

class linearInterpolation
:
    public topOInterpolationFunction
{

private:

    // Private Member Functions

        //- No copy construct
        linearInterpolation(const linearInterpolation&);

        //- No copy assignment
        void operator=(const linearInterpolation&);


public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from components
        linearInterpolation
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~linearInterpolation() = default;


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
