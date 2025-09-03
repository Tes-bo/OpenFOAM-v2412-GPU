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
    Foam::DiagonalMatrix

Description
    A templated (N x N) diagonal matrix of objects of \<Type\>, effectively
    containing N elements, derived from List.

See also
    Test-DiagonalMatrix.C

SourceFiles
    DiagonalMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_DiagonalMatrix_H
#define Foam_DiagonalMatrix_H

#include "List.H"
#include <numeric>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Form, class Type> class Matrix;

/*---------------------------------------------------------------------------*\
                       Class DiagonalMatrix Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class DiagonalMatrix
:
    public List<Type>
{
public:

    // Generated Methods

        //- Default construct
        DiagonalMatrix() = default;

        //- Copy construct
        DiagonalMatrix(const DiagonalMatrix&) = default;

        //- Copy assignment
        DiagonalMatrix& operator=(const DiagonalMatrix&) = default;


    // Constructors

        //- Construct from size, uninitialised content
        explicit DiagonalMatrix(const label n);

        //- Construct from size and initialise all elems to zero
        DiagonalMatrix(const label n, const Foam::zero);

        //- Construct from size and initialise all elems to value
        DiagonalMatrix(const label n, const Type& val);

        //- Construct from the diagonal of a Matrix
        template<class Form>
        DiagonalMatrix(const Matrix<Form, Type>& mat);


    // Member Functions

        //- Return the matrix inverse into itself
        void invert();

        //- Return a sort permutation labelList according to
        //- a given comparison on the diagonal entries
        template<class CompOp>
        List<label> sortPermutation(CompOp& compare) const;

        //- Column-reorder this Matrix according to
        //- a given permutation labelList
        void applyPermutation(const List<label>& p);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "DiagonalMatrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
