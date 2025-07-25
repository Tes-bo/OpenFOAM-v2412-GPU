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
    Global functions and operators related to the MinMax class.
    Included by MinMax.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_MinMaxOps_H
#define Foam_MinMaxOps_H

#include "MinMax.H"
#include "VectorSpace.H"
#include <type_traits>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Global Functions

//- The mag() function for min/max range
template<class T>
inline scalar mag(const MinMax<T>& range)
{
    return range.mag();
}


//- Return the value after clamping by the min/max limiter
template<class T>
inline T clip(const T& val, const MinMax<T>& range)
{
    return range.clamp(val);
}


//- Unary function for applying component-wise clamping
template<class T>
struct clampOp
{
    const T lower;
    const T upper;

    //- Construct from min/max limits. No validity checks
    clampOp(const T& min, const T& max)
    :
        lower(min),
        upper(max)
    {}

    //- Construct from min/max range. No validity checks
    clampOp(const MinMax<T>& range)
    :
        lower(range.min()),
        upper(range.max())
    {}

    //- Construct as 0-1 min/max range
    clampOp(const Foam::zero_one)
    :
        clampOp(MinMax<T>(Foam::zero_one{}))
    {}

    T operator()(const T& val) const
    {
        return min(max(val, lower), upper);
    }
};


/// Not really needed
/// //- Unary function for applying component-wise clamping (inplace)
/// template<class T>
/// struct clampEqOp : public clampOp<T>
/// {
///     using clampOp<T>::clampOp;
///
///     void operator()(T& val) const
///     {
///         val = clampOp<T>::operator()(val);
///     }
/// };


//- Extract the min/max range from a list of values.
template<class T>
inline MinMax<T> minMax(const UList<T>& vals)
{
    return MinMax<T>(vals);
}


//- Combine two values to create a min/max range. Order is unimportant.
template<class T>
inline MinMax<T> minMax(const T& x, const T& y)
{
    return MinMax<T>(x).add(y);
}


//- Combine two MinMax ranges (same as x + y)
template<class T>
inline MinMax<T> minMax(const MinMax<T>& x, const MinMax<T>& y)
{
    return MinMax<T>(x).add(y);
}


//- Combine values and/or MinMax ranges
template<class T>
struct minMaxOp
{
    MinMax<T> operator()(const T& x, const T& y) const
    {
        return MinMax<T>(x).add(y);
    }

    MinMax<T> operator()(const MinMax<T>& x, const T& y) const
    {
        return MinMax<T>(x).add(y);
    }

    MinMax<T> operator()(const T& x, const MinMax<T>& y) const
    {
        return MinMax<T>(y).add(x);
    }

    MinMax<T> operator()(const MinMax<T>& x, const MinMax<T>& y) const
    {
        return MinMax<T>(x).add(y);  // Same as (x + y)
    }
};


//- Combine assignment for MinMax range
template<class T>
struct minMaxEqOp
{
    MinMax<T>& operator()(MinMax<T>& x, const T& y) const
    {
        return x.add(y);
    }

    MinMax<T>& operator()(MinMax<T>& x, const UList<T>& y) const
    {
        return x.add(y);
    }

    MinMax<T>& operator()(MinMax<T>& x, const MinMax<T>& y) const
    {
        return x.add(y);
    }
};


//- The magnitude of a single value.
inline scalarMinMax minMaxMag(const scalar val)
{
    return scalarMinMax(Foam::mag(val));
}


//- The magnitude of a VectorSpace.
template<class Form, class Cmpt, direction nCmpt>
inline scalarMinMax minMaxMag(const VectorSpace<Form,Cmpt,nCmpt>& vs)
{
    return scalarMinMax(Foam::mag(vs));
}


//- The min/max magnitudes from a list of values
template<class T>
inline scalarMinMax minMaxMag(const UList<T>& vals)
{
    scalarMinMax result;
    for (const T& val : vals)
    {
        result += Foam::mag(val);
    }

    return result;
}


//- The min/max magnitudes from a min/max range
template<class T>
inline scalarMinMax minMaxMag(const MinMax<T>& range)
{
    return
    (
        scalarMinMax(Foam::mag(range.min())).add(Foam::mag(range.max()))
    );
}


//- Combine the magitude of two values to create a min/max range.
//- Order is unimportant.
template<class T>
inline scalarMinMax minMaxMag(const T& x, const T& y)
{
    return minMaxMag(x).add(Foam::mag(y));
}


//- Scalar combine two MinMax ranges of same type
template<class T>
inline scalarMinMax minMaxMag(const MinMax<T>& x, const MinMax<T>& y)
{
    return
    (
        minMaxMag(x)
        .add(Foam::mag(y.min()))
        .add(Foam::mag(y.max()))
    );
}


//- Scalar combine two MinMax ranges of dissimilar types
template<class T1, class T2>
inline scalarMinMax minMaxMag(const MinMax<T1>& x, const MinMax<T2>& y)
{
    return
    (
        minMaxMag(x)
        .add(Foam::mag(y.min()))
        .add(Foam::mag(y.max()))
    );
}


//- Scalar combine the magitude of a value.
template<class T>
struct minMaxMagOp
{
    scalarMinMax operator()(const scalarMinMax& x, const T& y) const
    {
        return minMaxMag(x).add(Foam::mag(y));
    }

    template<class T1, class T2>
    scalarMinMax operator()(const MinMax<T1>& x, const MinMax<T2>& y) const
    {
        return minMaxMag(x, y);
    }
};


//- Combine assignment for MinMax range
template<class T>
struct minMaxMagEqOp
{
    scalarMinMax& operator()(scalarMinMax& x, const T& y) const
    {
        x = minMaxMag(x);
        return x.add(Foam::mag(y));
    }

    scalarMinMax& operator()(scalarMinMax& x, const MinMax<T>& y) const
    {
        x = minMaxMag(x);

        return
        (
            x
            .add(Foam::mag(y.min()))
            .add(Foam::mag(y.max()))
        );
    }

    scalarMinMax& operator()(scalarMinMax& x, const UList<T>& y) const
    {
        x = minMaxMag(x);

        for (const T& val : y)
        {
            x.add(Foam::mag(val));
        }

        return x;
    }
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

//- Combine two ranges
template<class T>
inline MinMax<T> operator+(const MinMax<T>& x, const MinMax<T>& y)
{
    return MinMax<T>(x).add(y);
}


//- Multiply range by scalar factor
template<class T>
inline MinMax<T> operator*(const MinMax<T>& x, const scalar& s)
{
    return MinMax<T>(x.min()*s, x.max()*s);
}


//- Divide range by scalar factor
template<class T>
inline MinMax<T> operator/(const MinMax<T>& x, const scalar& s)
{
    return MinMax<T>(x.min()/s, x.max()/s);
}


// Comparison

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator<(const MinMax<T>& range, const U& val)
{
    return (range.compare(val) < 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator<=(const MinMax<T>& range, const U& val)
{
    return (range.compare(val) <= 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator>(const MinMax<T>& range, const U& val)
{
    return (range.compare(val) > 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator>=(const MinMax<T>& range, const U& val)
{
    return (range.compare(val) >= 0);
}


template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator<(const U& val, const MinMax<T>& range)
{
    return (range.compare(val) > 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator<=(const U& val, const MinMax<T>& range)
{
    return (range.compare(val) >= 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator>(const U& val, const MinMax<T>& range)
{
    return (range.compare(val) < 0);
}

template<class T, class U>
inline typename std::enable_if<std::is_convertible<U, T>::value, bool>::type
operator>=(const U& val, const MinMax<T>& range)
{
    return (range.compare(val) <= 0);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
