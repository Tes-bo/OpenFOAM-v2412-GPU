/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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
    Foam::vector

Description
    A Vector of values with scalar precision,
    where scalar is float/double depending on the compilation flags.

SourceFiles
    floatVectors.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_vector_H
#define Foam_vector_H

#include "scalar.H"
#include "Vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


//! \class Foam::floatVector
//! \brief A Vector of values with float precision
typedef Vector<float> floatVector;

//! \class Foam::doubleVector
//! \brief A Vector of values with double precision
typedef Vector<double> doubleVector;

// With float or double precision (depending on compilation)
typedef Vector<scalar> vector;

// With float or double precision (depending on compilation)
typedef Vector<solveScalar> solveVector;


// Traits

template<class Type>
class flux
:
    public innerProduct<vector, Type>
{};


template<>
class flux<scalar>
{
public:

    typedef scalar type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
