/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Namespace
    Foam::FieldOps

Description
    Various utility functions to work on Fields

SourceFiles
    FieldOps.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_FieldOps_H
#define Foam_FieldOps_H

#include "flipOp.H"
#include "ListOps.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace FieldOps
{

/*---------------------------------------------------------------------------*\
                        Namespace FieldOps Declarations
\*---------------------------------------------------------------------------*/

//- Populate a field as the result of a unary operation on an input.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
//
//  Examples,
//  \code
//    boolField result(sfield1.size());
//
//    FieldOps::assign(result, sfield1, lessOp1<scalar>(0.5));
//    FieldOps::assign(result, lfield1, std::logical_not<bool>());
//  \endcode
//
//  Example of using the Random::uniformGeneratorOp unary operator
//  to populate a random field,
//  \code
//    FieldOps::assign
//    (
//        sfield,
//        sfield,
//        Random::uniformGeneratorOp<scalar>(-15, 25)
//    );
//  \endcode
//
//  \note wraps std::transform
template<class Tout, class T1, class UnaryOp>
void assign
(
    Field<Tout>& result,
    const Field<T1>& a,
    const UnaryOp& op
);


//- Populate a field as the result of a binary operation on two inputs.
//  It is permissible for inputs/outputs to refer to the same field(s),
//  but partially overlapping regions are ill-defined.
//
//  Examples,
//  \code
//    FieldOps::assign(result, sfield1, sfield2, std::less<scalar>());
//    FieldOps::assign(result, lfield1, lfield2, std::logical_or<bool>());
//  \endcode
//
//  \note wraps std::transform
template<class Tout, class T1, class T2, class BinaryOp>
void assign
(
    Field<Tout>& result,
    const Field<T1>& a,
    const Field<T2>& b,
    const BinaryOp& bop
);


//- Emulate a ternary operation, selecting values from a or b
//- depending on the binary predicate.
//
//  Examples,
//  \code
//    FieldOps::ternary(result, sfield1, sfield2, std::less<scalar>());
//  \endcode
template<class T, class BinaryOp>
void ternary
(
    Field<T>& result,
    const Field<T>& a,
    const Field<T>& b,
    const BinaryOp& bop
);


//- Emulate a ternary operation, selecting values from a or b
//- depending on the conditional.
//
//  The meaning of the conditional is adjusted by the flip operation,
//  which is typically Foam::identityOp() or Foam::flipBoolOp()
//
//  Similar parameter requirements as Foam::subset()
//
//  Examples,
//  \code
//    FieldOps::ternarySelect(result, selector, sfield1, sfield2);
//  \endcode
template<class T, class BoolListType, class FlipOp>
void ternarySelect
(
    Field<T>& result,
    const BoolListType& cond,
    const Field<T>& a,
    const Field<T>& b,
    const FlipOp& flip
);


//- Emulate a ternary operation, selecting values from a or b
//- depending on the conditional.
//
//  The meaning of the conditional is adjusted by the flip operation,
//  which is typically Foam::identityOp() or Foam::flipBoolOp()
template<class T, class FlipOp>
void ternarySelect
(
    Field<T>& result,
    const bitSet& cond,
    const Field<T>& a,
    const Field<T>& b,
    const FlipOp& flip
);


//- Emulated ternary operation, without condition flipping
template<class T, class BoolListType>
void ternarySelect
(
    Field<T>& result,
    const BoolListType& cond,
    const Field<T>& a,
    const Field<T>& b
)
{
    ternarySelect(result, cond, a, b, identityOp{});
}


//- Emulated ternary operation, without condition flipping
template<class T>
void ternarySelect
(
    Field<T>& result,
    const bitSet& cond,
    const Field<T>& a,
    const Field<T>& b
)
{
    ternarySelect(result, cond, a, b, identityOp{});
}


//- Locate the min value in a field and return it and associated data
//  This can be used, for example, to return a value/position combination.
template<class T1, class T2>
Tuple2<T1,T2> findMinData
(
    const Field<T1>& vals,
    const Field<T2>& data
);


//- Locate the max value in a field and return it and associated data
//  This can be used, for example, to return a value/position combination.
template<class T1, class T2>
Tuple2<T1,T2> findMaxData
(
    const Field<T1>& vals,
    const Field<T2>& data
);


} // End namespace FieldOps

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "FieldOps.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
