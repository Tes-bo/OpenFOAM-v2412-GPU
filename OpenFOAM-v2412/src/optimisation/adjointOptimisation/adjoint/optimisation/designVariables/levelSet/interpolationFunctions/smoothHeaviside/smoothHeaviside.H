/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2007-2019 PCOpt/NTUA
    Copyright (C) 2013-2019 FOSS GP
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::smoothHeaviside

Description
    A smooth Heaviside function to project the signed distance field in 
    level set topology optimization.

\*---------------------------------------------------------------------------*/

#ifndef smoothHeaviside_H
#define smoothHeaviside_H

#include "topOInterpolationFunction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class smoothHeaviside Declaration
\*---------------------------------------------------------------------------*/

class smoothHeaviside
:
    public topOInterpolationFunction
{
protected:

    // Protected Data Members

        autoPtr<Function1<scalar>> b_;


private:
    
    // Private Member Functions

        //- No copy construct
        smoothHeaviside(const smoothHeaviside&) = delete;

        //- No copy assignment
        void operator=(const smoothHeaviside&) = delete;


public:

    //- Runtime type information
    TypeName("smoothHeaviside");


    // Constructors

        //- Construct from mesh and dictionary
        smoothHeaviside
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Destructor
    virtual ~smoothHeaviside() = default;


    // Member Functions

        //- Interpolate argument to result
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
