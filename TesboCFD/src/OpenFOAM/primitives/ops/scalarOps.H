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

InNamespace
    Foam

Description
    Functors that are scalar-specific.

\*---------------------------------------------------------------------------*/

#ifndef scalarOps_H
#define scalarOps_H

#include "scalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hypot operation (scalar only)
template<class T>
struct hypotOp
{
    T operator()(const T& x, const T& y) const
    {
        return std::hypot(x, y);
    }
};


//- Scalar division with divide-by-zero protection
//  Uses stabilise, but could also handle as per modulo and return zero
template<class T, class T2=Foam::scalar>
struct scalarDivideOp
{
    T operator()(const T& x, const T2& y) const
    {
        return (x / stabilise(y, pTraits<T2>::vsmall));
    }
};


//- Floating point modulo operation with divide-by-zero protection
template<class T, class T2=Foam::scalar>
struct scalarModuloOp
{
    T operator()(const T& x, const T2& y) const
    {
        if (Foam::mag(y) < pTraits<T2>::vsmall)
        {
            return pTraits<T>::zero;
        }
        return std::fmod(x, y);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
