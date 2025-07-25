/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    3D tensor transformation operations.

\*---------------------------------------------------------------------------*/

#ifndef Foam_transform_H
#define Foam_transform_H

#include "tensor.H"
#include "mathematicalConstants.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Rotational transformation tensor from vector n1 to n2
inline tensor rotationTensor
(
    const vector& n1,
    const vector& n2
)
{
    const scalar s = n1 & n2;
    const vector n3 = n1 ^ n2;
    const scalar magSqrN3 = magSqr(n3);

    // n1 and n2 define a plane n3
    if (magSqrN3 > SMALL)
    {
        // Return rotational transformation tensor in the n3-plane
        return
            s*I
          + (1 - s)*sqr(n3)/magSqrN3
          + (n2*n1 - n1*n2);
    }
    // n1 and n2 are contradirectional
    else if (s < 0)
    {
        // Return mirror transformation tensor
        return I + 2*n1*n2;
    }
    // n1 and n2 are codirectional
    else
    {
        // Return null transformation tensor
        return I;
    }
}


//- Rotational transformation tensor about the x-axis by omega radians
inline tensor Rx(const scalar omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        1,  0,  0,
        0,  c,  s,
        0, -s,  c
    );
}


//- Rotational transformation tensor about the y-axis by omega radians
inline tensor Ry(const scalar omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        c,  0, -s,
        0,  1,  0,
        s,  0,  c
    );
}


//- Rotational transformation tensor about the z-axis by omega radians
inline tensor Rz(const scalar omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);
    return tensor
    (
        c,  s,  0,
       -s,  c,  0,
        0,  0,  1
    );
}


//- Rotational transformation tensor about axis a by omega radians
inline tensor Ra(const vector& a, const scalar omega)
{
    const scalar s = sin(omega);
    const scalar c = cos(omega);

    return tensor
    (
        sqr(a.x())*(1 - c)  + c,
        a.y()*a.x()*(1 - c) + a.z()*s,
        a.x()*a.z()*(1 - c) - a.y()*s,

        a.x()*a.y()*(1 - c) - a.z()*s,
        sqr(a.y())*(1 - c)  + c,
        a.y()*a.z()*(1 - c) + a.x()*s,

        a.x()*a.z()*(1 - c) + a.y()*s,
        a.y()*a.z()*(1 - c) - a.x()*s,
        sqr(a.z())*(1 - c)  + c
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- No-op rotational transform for base types
template<class T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
transform(const tensor&, const T val)
{
    return val;
}


//- No-op rotational transform for spherical tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> transform
(
    const tensor&,
    const SphericalTensor<Cmpt>& val
)
{
    return val;
}


//- No-op inverse rotational transform for base types
template<class T>
constexpr typename std::enable_if<std::is_arithmetic<T>::value, T>::type
invTransform(const tensor&, const T val)
{
    return val;
}


//- No-op inverse rotational transform for spherical tensor
template<class Cmpt>
inline SphericalTensor<Cmpt> invTransform
(
    const tensor&,
    const SphericalTensor<Cmpt>& val
)
{
    return val;
}


//- Use rotational tensor to transform a vector.
//  Same as (rot & v)
template<class Cmpt>
inline Vector<Cmpt> transform(const tensor& tt, const Vector<Cmpt>& v)
{
    return (tt & v);
}


//- Use rotational tensor to inverse transform a vector.
//  Same as (v & rot)
template<class Cmpt>
inline Vector<Cmpt> invTransform(const tensor& tt, const Vector<Cmpt>& v)
{
    return (v & tt);
}


//- Use rotational tensor to transform a tensor.
//  Same as (rot & input & rot.T())
template<class Cmpt>
inline Tensor<Cmpt> transform(const tensor& tt, const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        // xx:
        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.xx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.xy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.xz(),

        // xy:
        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.yx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.yy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.yz(),

        // xz:
        (tt.xx()*t.xx() + tt.xy()*t.yx() + tt.xz()*t.zx())*tt.zx()
      + (tt.xx()*t.xy() + tt.xy()*t.yy() + tt.xz()*t.zy())*tt.zy()
      + (tt.xx()*t.xz() + tt.xy()*t.yz() + tt.xz()*t.zz())*tt.zz(),

        // yx:
        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.xx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.xy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.xz(),

        // yy:
        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.yx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.yy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.yz(),

        // yz:
        (tt.yx()*t.xx() + tt.yy()*t.yx() + tt.yz()*t.zx())*tt.zx()
      + (tt.yx()*t.xy() + tt.yy()*t.yy() + tt.yz()*t.zy())*tt.zy()
      + (tt.yx()*t.xz() + tt.yy()*t.yz() + tt.yz()*t.zz())*tt.zz(),

        // zx:
        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.xx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.xy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.xz(),

        // zy:
        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.yx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.yy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.yz(),

        // zz:
        (tt.zx()*t.xx() + tt.zy()*t.yx() + tt.zz()*t.zx())*tt.zx()
      + (tt.zx()*t.xy() + tt.zy()*t.yy() + tt.zz()*t.zy())*tt.zy()
      + (tt.zx()*t.xz() + tt.zy()*t.yz() + tt.zz()*t.zz())*tt.zz()
    );
}


//- Use rotational tensor to inverse transform a tensor.
//  Same as (rot.T() & input & rot)
template<class Cmpt>
inline Tensor<Cmpt> invTransform(const tensor& tt, const Tensor<Cmpt>& t)
{
    return Tensor<Cmpt>
    (
        // xx:
        (tt.xx()*t.xx() + tt.yx()*t.yx() + tt.zx()*t.zx())*tt.xx()
      + (tt.xx()*t.xy() + tt.yx()*t.yy() + tt.zx()*t.zy())*tt.yx()
      + (tt.xx()*t.xz() + tt.yx()*t.yz() + tt.zx()*t.zz())*tt.zx(),

        // xy:
        (tt.xx()*t.xx() + tt.yx()*t.yx() + tt.zx()*t.zx())*tt.xy()
      + (tt.xx()*t.xy() + tt.yx()*t.yy() + tt.zx()*t.zy())*tt.yy()
      + (tt.xx()*t.xz() + tt.yx()*t.yz() + tt.zx()*t.zz())*tt.zy(),

        // xz:
        (tt.xx()*t.xx() + tt.yx()*t.yx() + tt.zx()*t.zx())*tt.xz()
      + (tt.xx()*t.xy() + tt.yx()*t.yy() + tt.zx()*t.zy())*tt.yz()
      + (tt.xx()*t.xz() + tt.yx()*t.yz() + tt.zx()*t.zz())*tt.zz(),

        // yx:
        (tt.xy()*t.xx() + tt.yy()*t.yx() + tt.zy()*t.zx())*tt.xx()
      + (tt.xy()*t.xy() + tt.yy()*t.yy() + tt.zy()*t.zy())*tt.yx()
      + (tt.xy()*t.xz() + tt.yy()*t.yz() + tt.zy()*t.zz())*tt.zx(),

        // yy:
        (tt.xy()*t.xx() + tt.yy()*t.yx() + tt.zy()*t.zx())*tt.xy()
      + (tt.xy()*t.xy() + tt.yy()*t.yy() + tt.zy()*t.zy())*tt.yy()
      + (tt.xy()*t.xz() + tt.yy()*t.yz() + tt.zy()*t.zz())*tt.zy(),

        // yz:
        (tt.xy()*t.xx() + tt.yy()*t.yx() + tt.zy()*t.zx())*tt.xz()
      + (tt.xy()*t.xy() + tt.yy()*t.yy() + tt.zy()*t.zy())*tt.yz()
      + (tt.xy()*t.xz() + tt.yy()*t.yz() + tt.zy()*t.zz())*tt.zz(),

        // zx:
        (tt.xz()*t.xx() + tt.yz()*t.yx() + tt.zz()*t.zx())*tt.xx()
      + (tt.xz()*t.xy() + tt.yz()*t.yy() + tt.zz()*t.zy())*tt.yx()
      + (tt.xz()*t.xz() + tt.yz()*t.yz() + tt.zz()*t.zz())*tt.zx(),

        // zy:
        (tt.xz()*t.xx() + tt.yz()*t.yx() + tt.zz()*t.zx())*tt.xy()
      + (tt.xz()*t.xy() + tt.yz()*t.yy() + tt.zz()*t.zy())*tt.yy()
      + (tt.xz()*t.xz() + tt.yz()*t.yz() + tt.zz()*t.zz())*tt.zy(),

        // zz:
        (tt.xz()*t.xx() + tt.yz()*t.yx() + tt.zz()*t.zx())*tt.xz()
      + (tt.xz()*t.xy() + tt.yz()*t.yy() + tt.zz()*t.zy())*tt.yz()
      + (tt.xz()*t.xz() + tt.yz()*t.yz() + tt.zz()*t.zz())*tt.zz()
    );
}


//- Use rotational tensor to transform a symmetrical tensor.
//  Same as (rot & input & rot.T())
template<class Cmpt>
inline SymmTensor<Cmpt> transform(const tensor& tt, const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        // xx:
        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.xx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.xy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.xz(),

        // xy:
        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.yx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.yy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.yz(),

        // xz:
        (tt.xx()*st.xx() + tt.xy()*st.xy() + tt.xz()*st.xz())*tt.zx()
      + (tt.xx()*st.xy() + tt.xy()*st.yy() + tt.xz()*st.yz())*tt.zy()
      + (tt.xx()*st.xz() + tt.xy()*st.yz() + tt.xz()*st.zz())*tt.zz(),

        // yy:
        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.yx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.yy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.yz(),

        // yz:
        (tt.yx()*st.xx() + tt.yy()*st.xy() + tt.yz()*st.xz())*tt.zx()
      + (tt.yx()*st.xy() + tt.yy()*st.yy() + tt.yz()*st.yz())*tt.zy()
      + (tt.yx()*st.xz() + tt.yy()*st.yz() + tt.yz()*st.zz())*tt.zz(),

        // zz:
        (tt.zx()*st.xx() + tt.zy()*st.xy() + tt.zz()*st.xz())*tt.zx()
      + (tt.zx()*st.xy() + tt.zy()*st.yy() + tt.zz()*st.yz())*tt.zy()
      + (tt.zx()*st.xz() + tt.zy()*st.yz() + tt.zz()*st.zz())*tt.zz()
    );
}


//- Use rotational tensor to inverse transform a symmetrical tensor.
//  Same as (rot.T() & input & rot)
template<class Cmpt>
inline SymmTensor<Cmpt>
invTransform(const tensor& tt, const SymmTensor<Cmpt>& st)
{
    return SymmTensor<Cmpt>
    (
        // xx:
        (tt.xx()*st.xx() + tt.yx()*st.xy() + tt.zx()*st.xz())*tt.xx()
      + (tt.xx()*st.xy() + tt.yx()*st.yy() + tt.zx()*st.yz())*tt.yx()
      + (tt.xx()*st.xz() + tt.yx()*st.yz() + tt.zx()*st.zz())*tt.zx(),

        // xy:
        (tt.xx()*st.xx() + tt.yx()*st.xy() + tt.zx()*st.xz())*tt.xy()
      + (tt.xx()*st.xy() + tt.yx()*st.yy() + tt.zx()*st.yz())*tt.yy()
      + (tt.xx()*st.xz() + tt.yx()*st.yz() + tt.zx()*st.zz())*tt.zy(),

        // xz:
        (tt.xx()*st.xx() + tt.yx()*st.xy() + tt.zx()*st.xz())*tt.xz()
      + (tt.xx()*st.xy() + tt.yx()*st.yy() + tt.zx()*st.yz())*tt.yz()
      + (tt.xx()*st.xz() + tt.yx()*st.yz() + tt.zx()*st.zz())*tt.zz(),

        // yy:
        (tt.xy()*st.xx() + tt.yy()*st.xy() + tt.zy()*st.xz())*tt.xy()
      + (tt.xy()*st.xy() + tt.yy()*st.yy() + tt.zy()*st.yz())*tt.yy()
      + (tt.xy()*st.xz() + tt.yy()*st.yz() + tt.zy()*st.zz())*tt.zy(),

        // yz:
        (tt.xy()*st.xx() + tt.yy()*st.xy() + tt.zy()*st.xz())*tt.xz()
      + (tt.xy()*st.xy() + tt.yy()*st.yy() + tt.zy()*st.yz())*tt.yz()
      + (tt.xy()*st.xz() + tt.yy()*st.yz() + tt.zy()*st.zz())*tt.zz(),

        // zz:
        (tt.xz()*st.xx() + tt.yz()*st.xy() + tt.zz()*st.xz())*tt.xz()
      + (tt.xz()*st.xy() + tt.yz()*st.yy() + tt.zz()*st.yz())*tt.yz()
      + (tt.xz()*st.xz() + tt.yz()*st.yz() + tt.zz()*st.zz())*tt.zz()
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type1, class Type2>
inline Type1 transformMask(const Type2& t)
{
    return t;
}


template<>
inline sphericalTensor transformMask<sphericalTensor>(const tensor& t)
{
    return sph(t);
}


template<>
inline symmTensor transformMask<symmTensor>(const tensor& t)
{
    return symm(t);
}


//- Estimate angle of vec in coordinate system (e0, e1, e0^e1).
//  Is guaranteed to return increasing number but is not correct
//  angle. Used for sorting angles.  All input vectors need to be normalized.
//
//  Calculates scalar which increases with angle going from e0 to vec in
//  the coordinate system e0, e1, e0^e1
//
//  Jumps from 2*pi -> 0 at -SMALL so hopefully parallel vectors with small
//  rounding errors should still get the same quadrant.
//
inline scalar pseudoAngle
(
    const vector& e0,
    const vector& e1,
    const vector& vec
)
{
    const scalar cos_angle = vec & e0;
    const scalar sin_angle = vec & e1;

    if (sin_angle < -SMALL)
    {
        return (3.0 + cos_angle)*constant::mathematical::piByTwo;
    }
    else
    {
        return (1.0 - cos_angle)*constant::mathematical::piByTwo;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
