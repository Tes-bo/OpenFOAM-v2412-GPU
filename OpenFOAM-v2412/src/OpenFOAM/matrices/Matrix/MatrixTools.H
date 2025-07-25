/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Namespace
    Foam::MatrixTools

Description
    Collection of functions for matrix-related verifications.

SourceFiles
    MatrixTools.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_MatrixTools_H
#define Foam_MatrixTools_H

#include "Matrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Namespace MatrixTools Declaration
\*---------------------------------------------------------------------------*/

namespace MatrixTools
{

//- Compare matrix elements for absolute or relative equality
//
template<class Form1, class Form2, class Type>
bool equal
(
    const Matrix<Form1, Type>& A,
    const Matrix<Form2, Type>& B,
    const bool verbose = false,
    const label maxDiffs = 10, //!< Stop reporting after maxDiffs (0 to disable)
    const scalar relTol = 1e-5,
    const scalar absTol = 1e-8
);


//- Simple ASCII output of Matrix, MatrixBlock
template<class Container>
Ostream& printMatrix(Ostream& os, const Container& mat);


} // End namespace MatrixTools


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MatrixTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
