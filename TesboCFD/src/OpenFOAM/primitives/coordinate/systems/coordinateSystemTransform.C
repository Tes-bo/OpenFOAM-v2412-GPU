/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "coordinateSystem.H"
#include "transform.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Resolve templated global functions via local non-templated function.
// Lambda functions in the caller is a much messier solution.

#undef  makeTransform
#define makeTransform(Op, Type)                                               \
    static inline Type Op##_##Type(const tensor& tt, const Type& in)          \
    {                                                                         \
        return Op(tt, in);                                                    \
    }

makeTransform(transform, scalar);
makeTransform(transform, vector);
makeTransform(transform, sphericalTensor);
makeTransform(transform, symmTensor);
makeTransform(transform, tensor);

makeTransform(invTransform, scalar);
makeTransform(invTransform, vector);
makeTransform(invTransform, sphericalTensor);
makeTransform(invTransform, symmTensor);
makeTransform(invTransform, tensor);

#undef makeTransform

    //- Transform principal.
    static inline symmTensor transformPrincipal_vector
    (
        const tensor& tt,
        const vector& v
    )
    {
        return symmTensor
        (
            tt.xx()*v.x()*tt.xx()
          + tt.xy()*v.y()*tt.xy()
          + tt.xz()*v.z()*tt.xz(),

            tt.xx()*v.x()*tt.yx()
          + tt.xy()*v.y()*tt.yy()
          + tt.xz()*v.z()*tt.yz(),

            tt.xx()*v.x()*tt.zx()
          + tt.xy()*v.y()*tt.zy()
          + tt.xz()*v.z()*tt.zz(),

            tt.yx()*v.x()*tt.yx()
          + tt.yy()*v.y()*tt.yy()
          + tt.yz()*v.z()*tt.yz(),

            tt.yx()*v.x()*tt.zx()
          + tt.yy()*v.y()*tt.zy()
          + tt.yz()*v.z()*tt.zz(),

            tt.zx()*v.x()*tt.zx()
          + tt.zy()*v.y()*tt.zy()
          + tt.zz()*v.z()*tt.zz()
        );
    }

} // End namespace Foam


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::pointField> Foam::coordinateSystem::transformPoint
(
    const UList<point>& localCart
) const
{
    return transformPointImpl(localCart);
}


Foam::tmp<Foam::pointField> Foam::coordinateSystem::transformPoint
(
    const pointUIndList& localCart
) const
{
    return transformPointImpl(localCart);
}


Foam::tmp<Foam::pointField> Foam::coordinateSystem::invTransformPoint
(
    const UList<point>& global
) const
{
    return invTransformPointImpl(global);
}


Foam::tmp<Foam::pointField> Foam::coordinateSystem::invTransformPoint
(
    const pointUIndList& global
) const
{
    return invTransformPointImpl(global);
}


// Transformations

#undef  makeCoordinateSystemTransform
#define makeCoordinateSystemTransform(Op, RetType, Type)                      \
    Foam::RetType Foam::coordinateSystem::Op                                  \
    (                                                                         \
        const Type& input                                                     \
    ) const                                                                   \
    {                                                                         \
        return Op##_##Type(rot_, input);                                      \
    }                                                                         \
                                                                              \
    Foam::tmp<Foam::Field<Foam::RetType>> Foam::coordinateSystem::Op          \
    (                                                                         \
        const UList<Type>& input                                              \
    ) const                                                                   \
    {                                                                         \
        return manyTimesImpl<RetType>(rot_, input, Op##_##Type);              \
    }                                                                         \
                                                                              \
    Foam::RetType Foam::coordinateSystem::Op                                  \
    (                                                                         \
        const point& global,                                                  \
        const Type& input                                                     \
    ) const                                                                   \
    {                                                                         \
        return Op##_##Type(this->R(global), input);                           \
    }                                                                         \
                                                                              \
    Foam::tmp<Foam::Field<Foam::RetType>> Foam::coordinateSystem::Op          \
    (                                                                         \
        const UList<point>& global,                                           \
        const Type& input                                                     \
    ) const                                                                   \
    {                                                                         \
        return oneToManyImpl<RetType>(global, input, Op##_##Type);            \
    }                                                                         \
                                                                              \
    Foam::tmp<Foam::Field<Foam::RetType>> Foam::coordinateSystem::Op          \
    (                                                                         \
        const pointUIndList& global,                                          \
        const Type& input                                                     \
    ) const                                                                   \
    {                                                                         \
        return oneToManyImpl<RetType>(global, input, Op##_##Type);            \
    }                                                                         \
                                                                              \
    Foam::tmp<Foam::Field<Foam::RetType>> Foam::coordinateSystem::Op          \
    (                                                                         \
        const UList<point>& global,                                           \
        const UList<Type>& input                                              \
    ) const                                                                   \
    {                                                                         \
        return oneToOneImpl<RetType>(global, input, Op##_##Type);             \
    }                                                                         \
                                                                              \
    Foam::tmp<Foam::Field<Foam::RetType>> Foam::coordinateSystem::Op          \
    (                                                                         \
        const pointUIndList& global,                                          \
        const UList<Type>& input                                              \
    ) const                                                                   \
    {                                                                         \
        return oneToOneImpl<RetType>(global, input, Op##_##Type);             \
    }


makeCoordinateSystemTransform(transformPrincipal, symmTensor, vector);

makeCoordinateSystemTransform(transform, scalar, scalar);
makeCoordinateSystemTransform(transform, vector, vector);
makeCoordinateSystemTransform(transform, sphericalTensor, sphericalTensor);
makeCoordinateSystemTransform(transform, symmTensor, symmTensor);
makeCoordinateSystemTransform(transform, tensor, tensor);

makeCoordinateSystemTransform(invTransform, scalar, scalar);
makeCoordinateSystemTransform(invTransform, vector, vector);
makeCoordinateSystemTransform(invTransform, sphericalTensor, sphericalTensor);
makeCoordinateSystemTransform(invTransform, symmTensor, symmTensor);
makeCoordinateSystemTransform(invTransform, tensor, tensor);

#undef makeCoordinateSystemTransform


// ************************************************************************* //
