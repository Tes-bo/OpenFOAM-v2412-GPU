/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Typedef
    Foam::Scalar

Description
    Floating-point number (float or double)

Note
    The floor/ceil/round operations could easily be extended to handle
    VectorSpace when the need arises. For example,

    \code
    template<class T>
    struct floorOp
    {
        T operator()(const T& x) const
        {
            T ret;
            for (direction cmpt=0; cmpt < pTraits<T>::nComponents; ++cmpt)
            {
                setComponent(ret, cmpt) = std::floor(component(x, cmpt));
            }
            return ret;
        }
    };
    \endcode

SourceFiles
    Scalar.C

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Template specialisation for pTraits<Scalar>
template<>
class pTraits<Scalar>
{
    Scalar p_;

public:

    // Typedefs

        //- Component type
        typedef Scalar cmptType;

        //- Magnitude type
        typedef Scalar magType;

        //- Equivalent type of labels used for valid component indexing
        typedef label labelType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of Scalar is 0
        static constexpr direction rank = 0;

        //- Number of components in Scalar is 1
        static constexpr direction nComponents = 1;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];
        static const Scalar zero;
        static const Scalar one;
        static const Scalar max;
        static const Scalar min;
        static const Scalar rootMax;
        static const Scalar rootMin;
        static const Scalar vsmall;


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(Scalar val) noexcept
        :
            p_(val)
        {}

        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator Scalar() const noexcept
        {
            return p_;
        }

        //- Access the value
        operator Scalar&() noexcept
        {
            return p_;
        }
};


// * * * * * * * * * * * * * * * IO/Conversion * * * * * * * * * * * * * * * //

//- A word representation of a floating-point value.
//  Uses stringstream instead of std::to_string for more consistent formatting.
word name(const Scalar val);


//- A word representation of a floating-point value.
template<>
struct nameOp<Scalar>
{
    word operator()(const Scalar val) const
    {
        return Foam::name(val);
    }
};


//- Parse entire buffer as a float/double, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
Scalar ScalarRead(const char* buf);

//- Parse entire buffer as a float/double, skipping leading/trailing whitespace.
//  \return True if successful.
bool ScalarRead(const char* buf, Scalar& val);

//- Parse entire string as a float/double, skipping leading/trailing whitespace.
//  \return Parsed value or FatalIOError on any problem
inline Scalar ScalarRead(const std::string& str)
{
    return ScalarRead(str.c_str());
}

//- Parse entire string as a float/double, skipping leading/trailing whitespace.
//  \return True if successful.
inline bool ScalarRead(const std::string& str, Scalar& val)
{
    return ScalarRead(str.c_str(), val);
}


Scalar ScalarRead(Istream& is);
Istream& operator>>(Istream& is, Scalar& val);
Ostream& operator<<(Ostream& os, const Scalar val);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Standard C++ transcendental functions
transFunc(sqrt)

transFunc(cbrt)
transFunc(exp)
transFunc(log)
transFunc(log10)
transFunc(sin)
transFunc(cos)
transFunc(tan)
transFunc(asin)
transFunc(acos)
transFunc(atan)
transFunc(sinh)
transFunc(cosh)
transFunc(tanh)
transFunc(asinh)
transFunc(acosh)
transFunc(atanh)

// Standard ANSI-C (but not in <cmath>) transcendental functions

transFunc(erf)
transFunc(erfc)
transFunc(lgamma)
transFunc(tgamma)

besselFunc(j0)
besselFunc(j1)
besselFunc(y0)
besselFunc(y1)
besselFunc2(jn)
besselFunc2(yn)


//- Non-const access to scalar-type (has no components)
inline Scalar& setComponent(Scalar& val, const direction) noexcept
{
    return val;
}


//- Return scalar value (has no components)
inline constexpr Scalar component(const Scalar val, const direction) noexcept
{
    return val;
}


//- Return 1 if s is greater_equal zero, or otherwise -1
inline Scalar sign(const Scalar s) noexcept
{
    return (s >= 0)? 1: -1;
}


//- Return 1 if s is greater than zero, otherwise 0
inline Scalar pos(const Scalar s) noexcept
{
    return (s > 0)? 1: 0;
}


//- Return 1 if s is greater_equal zero, or otherwise 0
inline Scalar pos0(const Scalar s) noexcept
{
    return (s >= 0)? 1: 0;
}


//- Return 1 if s is less than zero, or otherwise 0
inline Scalar neg(const Scalar s) noexcept
{
    return (s < 0)? 1: 0;
}


//- Return 1 if s is less_equal zero, or otherwise 0
inline Scalar neg0(const Scalar s) noexcept
{
    return (s <= 0)? 1: 0;
}


//- Return the positive part of s, otherwise zero. Same as max(0, s).
inline Scalar posPart(const Scalar s) noexcept
{
    return (s > 0)? s: 0;
}


//- Return the negative part of s, otherwise zero. Same as min(0, s).
//  Does not change the sign
inline Scalar negPart(const Scalar s) noexcept
{
    return (s < 0)? s: 0;
}


inline bool equal(const Scalar& a, const Scalar& b)
{
    return mag(a - b) <= ScalarVSMALL;
}


inline bool notEqual(const Scalar a, const Scalar b)
{
    return mag(a - b) > ScalarVSMALL;
}


//- Clamp scalar value to a 0-1 range
inline Scalar clamp(const Scalar val, const Foam::zero_one)
{
    return (val < 0) ? 0 : (1 < val) ? 1 : val;
}


// Fast implementation, and with scalar promotion of upper/lower limits
inline Scalar clamp(const Scalar val, const Scalar lower, const Scalar upper)
{
    return (val < lower) ? lower : (upper < val) ? upper : val;
}


inline Scalar limit(const Scalar s1, const Scalar s2)
{
    return (mag(s1) < mag(s2)) ? s1: 0.0;
}


inline Scalar minMod(const Scalar s1, const Scalar s2)
{
    return (mag(s1) < mag(s2)) ? s1: s2;
}


//- Linear interpolation of scalars a and b by factor t
//  Explicit form (not fused multiply add etc) for exact values at end points.
//  Unlike std::lerp (C++20) there are no isfinite() checks etc.
inline constexpr Scalar lerp(const Scalar a, const Scalar b, const Scalar t)
{
    return (Scalar{1}-t)*a + t*b;
}


inline Scalar magSqr(const Scalar s)
{
    return s*s;
}


inline Scalar sqr(const Scalar s)
{
    return s*s;
}


inline Scalar pow3(const Scalar s)
{
    return s*sqr(s);
}


inline Scalar pow4(const Scalar s)
{
    return sqr(sqr(s));
}


inline Scalar pow5(const Scalar s)
{
    return s*pow4(s);
}


inline Scalar pow6(const Scalar s)
{
    return pow3(sqr(s));
}


inline Scalar pow025(const Scalar s)
{
    return sqrt(sqrt(s));
}


inline Scalar inv(const Scalar s)
{
    return 1.0/s;
}


inline Scalar dot(const Scalar s1, const Scalar s2)
{
    return s1*s2;
}


inline Scalar cmptMultiply(const Scalar s1, const Scalar s2)
{
    return s1*s2;
}


inline Scalar cmptPow(const Scalar s1, const Scalar s2)
{
    return pow(s1, s2);
}


inline Scalar cmptDivide(const Scalar s1, const Scalar s2)
{
    return s1/s2;
}


inline Scalar cmptMax(const Scalar s)
{
    return s;
}


inline Scalar cmptMin(const Scalar s)
{
    return s;
}


inline Scalar cmptAv(const Scalar s)
{
    return s;
}


inline Scalar cmptSqr(const Scalar s)
{
    return sqr(s);
}


inline Scalar cmptMag(const Scalar s)
{
    return mag(s);
}


inline Scalar cmptMagSqr(const Scalar s)
{
    return magSqr(s);
}


inline Scalar sqrtSumSqr(const Scalar a, const Scalar b)
{
    const Scalar maga = mag(a);
    const Scalar magb = mag(b);

    if (maga > magb)
    {
        return maga*sqrt(1.0 + sqr(magb/maga));
    }
    else
    {
        return magb < ScalarVSMALL ? 0.0 : magb*sqrt(1.0 + sqr(maga/magb));
    }
}


//- Stabilisation around zero for division
inline Scalar stabilise(const Scalar s, const Scalar tol)
{
    if (s >= 0)
    {
        return s + tol;
    }
    else
    {
        return s - tol;
    }
}


// Specializations

// Default definition in ops.H
template<class T> struct compareOp;

//- Compare scalar values
template<>
struct compareOp<Scalar>
{
    const Scalar tolerance;

    //- Construct with specified tolerance (non-negative value)
    compareOp(Scalar tol = ScalarVSMALL)
    :
        tolerance(tol)
    {}

    Scalar operator()(const Scalar& a, const Scalar& b) const
    {
        return (mag(a - b) <= tolerance) ? 0 : (a - b);
    }
};


// Default definition in ops.H
template<class T> struct equalOp;

//- Compare scalar values for equality
template<>
struct equalOp<Scalar>
{
    const Scalar tolerance;

    //- Construct with specified tolerance (non-negative value)
    equalOp(Scalar tol = ScalarVSMALL)
    :
        tolerance(tol)
    {}

    bool operator()(const Scalar& a, const Scalar& b) const
    {
        return mag(a - b) <= tolerance;
    }
};


// Default definition in ops.H
template<class T> struct notEqualOp;

//- Compare scalar values for inequality
template<>
struct notEqualOp<Scalar>
{
    const Scalar tolerance;

    //- Construct with specified tolerance (non-negative value)
    notEqualOp(Scalar tol = ScalarVSMALL)
    :
        tolerance(tol)
    {}

    bool operator()(const Scalar& a, const Scalar& b) const
    {
        return mag(a - b) > tolerance;
    }
};



// Default definition in ops.H (future?)
template<class T> struct floorOp;

//- Round scalar downwards - functor version of std::floor
template<>
struct floorOp<Scalar>
{
    Scalar operator()(const Scalar& x) const
    {
        return std::floor(x);
    }
};


// Default definition in ops.H (future?)
template<class T> struct ceilOp;

//- Round scalar upwards - functor version of std::ceil
template<>
struct ceilOp<Scalar>
{
    Scalar operator()(const Scalar& x) const
    {
        return std::ceil(x);
    }
};


// Default definition in ops.H (future?)
template<class T> struct roundOp;

//- Round scalar - functor version of std::round
template<>
struct roundOp<Scalar>
{
    Scalar operator()(const Scalar& x) const
    {
        return std::round(x);
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
