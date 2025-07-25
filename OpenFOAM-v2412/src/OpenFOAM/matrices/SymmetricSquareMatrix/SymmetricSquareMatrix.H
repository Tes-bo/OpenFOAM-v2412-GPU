/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::SymmetricSquareMatrix

Description
    A templated (N x N) square matrix of objects of \<Type\>,
    containing N*N elements, derived from Matrix.

See also
    Test-SymmetricSquareMatrix.C

SourceFiles
    SymmetricSquareMatrixI.H
    SymmetricSquareMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SymmetricSquareMatrix_H
#define Foam_SymmetricSquareMatrix_H

#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class SymmetricSquareMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class SymmetricSquareMatrix
:
    public Matrix<SymmetricSquareMatrix<Type>, Type>
{
public:

    // Generated Methods

        //- Default construct
        SymmetricSquareMatrix() = default;

        //- Copy construct
        SymmetricSquareMatrix(const SymmetricSquareMatrix&) = default;

        //- Copy assignment
        SymmetricSquareMatrix&
        operator=(const SymmetricSquareMatrix&) = default;


    // Constructors

        //- Construct for given size (rows == cols), uninitialised content
        inline explicit SymmetricSquareMatrix(const label n);

        //- Construct for given size (rows == cols)
        //- initializing all elements to zero
        inline SymmetricSquareMatrix(const label n, const Foam::zero);

        //- Construct for given size (rows == cols)
        //- initializing all elements to the given value
        inline SymmetricSquareMatrix(const label n, const Type& val);

        //- Construct for given size (rows == cols)
        //- initializing to the identity matrix
        template<class AnyType>
        inline SymmetricSquareMatrix(const label n, const Identity<AnyType>);

        //- Construct from Istream
        inline explicit SymmetricSquareMatrix(Istream& is);

        //- Clone
        inline autoPtr<SymmetricSquareMatrix<Type>> clone() const;


    // Member Operators

        //- Assign all elements to zero
        inline void operator=(const Foam::zero);

        //- Assign all elements to value
        inline void operator=(const Type& val);

        //- Set to identity matrix
        template<class AnyType>
        void operator=(const Identity<AnyType>);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Global Functions

//- Return the LU decomposed SymmetricSquareMatrix inverse
template<class Type>
SymmetricSquareMatrix<Type> invDecomposed(const SymmetricSquareMatrix<Type>&);

//- Return the SymmetricSquareMatrix inverse
template<class Type>
SymmetricSquareMatrix<Type> inv(const SymmetricSquareMatrix<Type>&);

//- Return the LU decomposed SymmetricSquareMatrix det
template<class Type>
Type detDecomposed(const SymmetricSquareMatrix<Type>&);

//- Return the SymmetricSquareMatrix det
template<class Type>
Type det(const SymmetricSquareMatrix<Type>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SymmetricSquareMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SymmetricSquareMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
