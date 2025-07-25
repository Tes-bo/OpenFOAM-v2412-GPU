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
    topOInterpolationFunction

Description
    Abstract base class for interpolation functions used in topology
    optimisation

SourceFiles
    topOInterpolationFunction.C

\*---------------------------------------------------------------------------*/

#ifndef topOInterpolationFunction_H
#define topOInterpolationFunction_H

#include "fvMesh.H"
#include "Function1.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class topOInterpolationFunction Declaration
\*---------------------------------------------------------------------------*/

class topOInterpolationFunction
{

protected:

    // Protected data

        const fvMesh& mesh_;

        dictionary dict_;


private:

    // Private Member Functions

        //- No copy construct
        topOInterpolationFunction(const topOInterpolationFunction&);

        //- No copy assignment
        void operator=(const topOInterpolationFunction&);


public:

    //- Runtime type information
    TypeName("topOInterpolationFunction");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            topOInterpolationFunction,
            dictionary,
            (
                const fvMesh& mesh,
                const dictionary& dict
            ),
            (mesh, dict)
        );


    // Constructors

        //- Construct from components
        topOInterpolationFunction
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    // Selectors

        //- Return an autoPtr to the selected interpolation type
        static autoPtr<topOInterpolationFunction> New
        (
            const fvMesh& mesh,
            const dictionary& dict
        );


    //- Destructor
    virtual ~topOInterpolationFunction() = default;


    // Member Functions

        //- Interpolate argument to result
        virtual void interpolate
        (
            const scalarField& arg,
            scalarField& res
        ) const = 0;

        //- Return of function with respect to the argument field
        virtual tmp<scalarField> derivative(const scalarField& arg) const = 0;

        //- Set the parameter determining length scale
        //  Warning if the selected interpolation function has no such param
        virtual void setLengthScaleParam(const scalar lengthScale);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
