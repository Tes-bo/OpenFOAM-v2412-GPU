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
    Foam::fa::leastSquaresFaGrad

Description
    Second-order gradient scheme using least-squares.

Author
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    leastSquaresFaGrad.C

\*---------------------------------------------------------------------------*/

#ifndef leastSquaresFaGrad_H
#define leastSquaresFaGrad_H

#include "faGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fa
{

/*---------------------------------------------------------------------------*\
                     Class leastSquaresFaGrad Declaration\*-------------------
\*---------------------------------------------------------------------------*/

template<class Type>
class leastSquaresFaGrad
:
    public fa::gradScheme<Type>
{
public:

    //- Runtime type information
    TypeName("leastSquares");


    // Generated Methods

        //- No copy construct
        leastSquaresFaGrad(const leastSquaresFaGrad&) = delete;

        //- No copy assignment
        void operator=(const leastSquaresFaGrad&) = delete;


    // Constructors

        //- Construct from mesh
        leastSquaresFaGrad(const faMesh& mesh)
        :
            gradScheme<Type>(mesh)
        {}

        //- Construct from Istream
        leastSquaresFaGrad(const faMesh& mesh, Istream&)
        :
            gradScheme<Type>(mesh)
        {}


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

#ifdef NoRepository
    #include "leastSquaresFaGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
