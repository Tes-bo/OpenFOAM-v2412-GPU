/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2017 Wikki Ltd
    Copyright (C) 2023 OpenCFD Ltd.
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
    Foam::fa::edgeLimitedGrad

Description
    edgeLimitedGrad gradient scheme applied to a runTime selected base gradient
    scheme.

    The scalar limiter based on limiting the extrapolated edge values
    between the edge-neighbour cell values and is applied to all components
    of the gradient.

SourceFiles
    edgeLimitedFaGrads.C

\*---------------------------------------------------------------------------*/

#ifndef edgeLimitedFaGrad_H
#define edgeLimitedFaGrad_H

#include "faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                       Class edgeLimitedGrad Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class edgeLimitedGrad
:
    public fa::gradScheme<Type>
{
    // Private Data

        //- Basic gradient scheme
        tmp<fa::gradScheme<Type>> basicGradScheme_;

        //- Limiter coefficient
        scalar k_;


    // Private Member Functions

        inline void limitEdge
        (
            scalar& limiter,
            const scalar maxDelta,
            const scalar minDelta,
            const scalar extrapolate
        ) const;


public:

    //- RunTime type information
    TypeName("edgeLimited");


    // Generated Methods

        //- No copy construct
        edgeLimitedGrad(const edgeLimitedGrad&) = delete;

        //- No copy assignment
        void operator=(const edgeLimitedGrad&) = delete;


    // Constructors

        //- Construct from mesh and schemeData
        edgeLimitedGrad(const faMesh& mesh, Istream& schemeData)
        :
            gradScheme<Type>(mesh),
            basicGradScheme_(fa::gradScheme<Type>::New(mesh, schemeData)),
            k_(readScalar(schemeData))
        {
            if (k_ < 0 || k_ > 1)
            {
                FatalIOErrorInFunction(schemeData)
                    << "coefficient = " << k_
                    << " should be >= 0 and <= 1"
                    << exit(FatalIOError);
            }
        }


    // Member Functions

        //- Return the gradient of the given field to the
        //- gradScheme::grad for optional caching
        virtual tmp
        <
            GeometricField
            <typename outerProduct<vector, Type>::type, faPatchField, areaMesh>
        > calcGrad
        (
            const GeometricField<Type, faPatchField, areaMesh>& vsf,
            const word& name
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fa

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
