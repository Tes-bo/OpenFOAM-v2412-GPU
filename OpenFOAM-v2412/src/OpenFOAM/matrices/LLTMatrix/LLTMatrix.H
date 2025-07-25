/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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
    Foam::LLTMatrix

Description
    Templated class to perform the Cholesky decomposition on a
    symmetric positive-definite matrix.

    Member functions are provided to solve linear systems using the LLT
    decomposition.

SourceFiles
    LLTMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef LLTMatrix_H
#define LLTMatrix_H

#include "SquareMatrix.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class LLTMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class LLTMatrix
:
    public SquareMatrix<Type>
{

    // Private Member Functions

        //- Solve the linear system with the given source
        //- and return the solution in x
        //  This function may be called with the same field for x and source.
        template<template<typename> class ListContainer>
        void solveImpl
        (
            List<Type>& x,
            const ListContainer<Type>& source
        ) const;


public:

    // Constructors

        //- Construct null
        LLTMatrix();

        //- Construct and perform the decomposition on input square matrix
        LLTMatrix(const SquareMatrix<Type>& mat);


    // Member Functions

        //- Copy matrix and perform Cholesky decomposition
        void decompose(const SquareMatrix<Type>& mat);

        //- Solve the linear system with the given source
        //- and return the solution in the argument x.
        //  This function may be called with the same field for x and source.
        void solve
        (
            List<Type>& x,
            const UList<Type>& source
        ) const;

        //- Solve the linear system with the given source
        //- and return the solution in the argument x.
        //  This function may be called with the same field for x and source.
        template<class Addr>
        void solve
        (
            List<Type>& x,
            const IndirectListBase<Type, Addr>& source
        ) const;

        //- Solve the linear system with the given source
        //- return the solution
        tmp<Field<Type>> solve
        (
            const UList<Type>& source
        ) const;

        //- Solve the linear system with the given source
        //- return the solution
        template<class Addr>
        tmp<Field<Type>> solve
        (
            const IndirectListBase<Type, Addr>& source
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "LLTMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
