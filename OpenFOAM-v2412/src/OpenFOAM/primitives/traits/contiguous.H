/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2019 OpenCFD Ltd.
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
    Foam::is_contiguous

Description
    A template class to specify that a data type can be considered as being
    contiguous in memory.

    Normally only integral and floating-point types can be considered
    contiguous, but some other types (eg, FixedList, Pair, Vector etc)
    consisting purely of these fundamental types can be considered
    as having a contiguous memory layout as well.

Note
    In OpenFOAM 1906 and earlier, the contiguous trait was handled
    by templated \c contiguous global functions.

    While possible to mark this as deleted, this does not detect or
    prevent specializations. Thus omit the usual housekeeping.

Class
    Foam::is_contiguous_label

Description
    A template class to specify if a data type is composed solely of
    Foam::label elements.

Class
    Foam::is_contiguous_scalar

Description
    A template class to specify if a data type is composed solely of
    Foam::scalar elements.

\*---------------------------------------------------------------------------*/

#ifndef Foam_contiguous_H
#define Foam_contiguous_H

#include "labelFwd.H"
#include "scalarFwd.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Base definition for (integral | floating-point) as contiguous
template<class T>
struct is_contiguous : std::is_arithmetic<T> {};


// Base definition for 'label'
template<class T>
struct is_contiguous_label : std::is_same<T, label> {};


// Base definition for 'scalar'
template<class T>
struct is_contiguous_scalar : std::is_same<T, scalar> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
