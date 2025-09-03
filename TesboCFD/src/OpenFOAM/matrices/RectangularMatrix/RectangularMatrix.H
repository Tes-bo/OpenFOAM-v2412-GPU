/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::RectangularMatrix

Description
    A templated (M x N) rectangular matrix of objects of \<Type\>,
    containing M*N elements, derived from Matrix.

See also
    Test-RectangularMatrix.C

SourceFiles
    RectangularMatrixI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_RectangularMatrix_H
#define Foam_RectangularMatrix_H

#include "Matrix.H"
#include "SquareMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class RectangularMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class RectangularMatrix
:
    public Matrix<RectangularMatrix<Type>, Type>
{
public:

    // Generated Methods

        //- Default construct
        RectangularMatrix() = default;

        //- Copy construct
        RectangularMatrix(const RectangularMatrix&) = default;

        //- Copy assignment
        RectangularMatrix& operator=(const RectangularMatrix&) = default;


    // Constructors

        //- Construct a square matrix (rows == columns), uninitialised content
        inline explicit RectangularMatrix(const label n);

        //- Construct given number of rows/columns
        inline RectangularMatrix(const label m, const label n);

        //- Construct given number of rows/columns
        //- initializing all elements to zero
        inline RectangularMatrix
        (
            const label m,
            const label n,
            const Foam::zero
        );

        //- Construct given number of rows/columns
        //- initializing all elements to the given value
        inline RectangularMatrix(const label m, const label n, const Type& val);

        //- Construct for given number of rows/columns
        //- initializing all elements to zero, and diagonal to one
        template<class AnyType>
        inline RectangularMatrix
        (
            const labelPair& dims,
            const Identity<AnyType>
        );

        //- Construct given number of rows/columns by using a label pair
        inline explicit RectangularMatrix(const labelPair& dims);

        //- Construct given number of rows/columns by using a label pair
        //- and initializing all elements to zero
        inline RectangularMatrix(const labelPair& dims, const Foam::zero);

        //- Construct given number of rows/columns by using a label pair
        //- and initializing all elements to the given value
        inline RectangularMatrix(const labelPair& dims, const Type& val);

        //- Construct from a block of another matrix
        template<class MatrixType>
        inline RectangularMatrix(const ConstMatrixBlock<MatrixType>& mat);

        //- Construct from a block of another matrix
        template<class MatrixType>
        inline RectangularMatrix(const MatrixBlock<MatrixType>& mat);

        //- Construct as copy of a square matrix
        inline RectangularMatrix(const SquareMatrix<Type>& mat);

        //- Construct from Istream
        inline explicit RectangularMatrix(Istream& is);

        //- Clone
        inline autoPtr<RectangularMatrix<Type>> clone() const;


    // Member Operators

        //- Move assignment
        inline void operator=(RectangularMatrix<Type>&& mat);

        //- Assign all elements to zero
        inline void operator=(const Foam::zero);

        //- Assign all elements to value
        inline void operator=(const Type& val);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "RectangularMatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
