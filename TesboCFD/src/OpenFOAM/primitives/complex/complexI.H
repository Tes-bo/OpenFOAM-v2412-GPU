/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2014 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline constexpr Foam::complex::complex() noexcept
:
    re(0),
    im(0)
{}


inline constexpr Foam::complex::complex(const Foam::zero) noexcept
:
    re(0),
    im(0)
{}


inline constexpr Foam::complex::complex(const scalar r) noexcept
:
    re(r),
    im(0)
{}


inline constexpr Foam::complex::complex(const scalar r, const scalar i) noexcept
:
    re(r),
    im(i)
{}


inline Foam::complex::complex(const std::complex<float>& c)
:
    re(c.real()),
    im(c.imag())
{}


inline Foam::complex::complex(const std::complex<double>& c)
:
    re(c.real()),
    im(c.imag())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::complex Foam::complex::conjugate() const
{
    return complex(re, -im);
}


inline Foam::scalar Foam::complex::magnitude() const
{
    return std::hypot(re, im);
}


inline Foam::scalar Foam::complex::magSqr() const
{
    return (re*re + im*im);
}


inline Foam::scalar Foam::complex::cmptSum() const noexcept
{
    return (re + im);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

inline void Foam::complex::operator=(const Foam::zero)
{
    re = 0;
    im = 0;
}


inline void Foam::complex::operator=(const scalar s)
{
    re = s;
    im = 0;
}


inline void Foam::complex::operator+=(const complex& c)
{
    re += c.re;
    im += c.im;
}


inline void Foam::complex::operator+=(const scalar s)
{
    re += s;
}


inline void Foam::complex::operator-=(const complex& c)
{
    re -= c.re;
    im -= c.im;
}


inline void Foam::complex::operator-=(const scalar s)
{
    re -= s;
}


inline void Foam::complex::operator*=(const complex& c)
{
    *this = (*this)*c;
}


inline void Foam::complex::operator*=(const scalar s)
{
    re *= s;
    im *= s;
}


inline void Foam::complex::operator/=(const complex& c)
{
    *this = *this/c;
}


inline void Foam::complex::operator/=(const scalar s)
{
    re /= s;
    im /= s;
}


inline bool Foam::complex::operator==(const complex& c) const
{
    return (equal(re, c.re) && equal(im, c.im));
}


inline bool Foam::complex::operator!=(const complex& c) const
{
    return !operator==(c);
}


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

inline Foam::complex Foam::operator~(const complex& c)
{
    return c.conjugate();
}


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

namespace Foam
{

inline scalar mag(const complex& c)
{
    return c.magnitude();
}


inline scalar magSqr(const complex& c)
{
    return c.magSqr();
}


inline complex sqr(const complex& c)
{
    return c*c;
}


//- sgn() https://en.wikipedia.org/wiki/Sign_function#Complex_signum
inline complex sign(const complex& c)
{
    // TBD: Use volatile to avoid aggressive branch optimization
    const scalar s(mag(c));
    return s < ROOTVSMALL ? Zero : c/s;
}


//- csgn() https://en.wikipedia.org/wiki/Sign_function#Complex_signum
inline scalar csign(const complex& c)
{
    return equal(c.real(), 0) ? sign(c.imag()) : sign(c.real());
}


inline const complex& min(const complex& c1, const complex& c2)
{
    return (c1.magSqr() < c2.magSqr()) ? c1 : c2;
}


inline const complex& max(const complex& c1, const complex& c2)
{
    return (c1.magSqr() < c2.magSqr()) ? c2 : c1;
}


inline complex limit(const complex& c1, const complex& c2)
{
    return complex
    (
        limit(c1.real(), c2.real()),
        limit(c1.imag(), c2.imag())
    );
}


//- Linear interpolation of complex a and b by factor t
inline complex lerp(const complex& a, const complex& b, const scalar t)
{
    const scalar onet = (1-t);

    return complex
    (
        onet*a.real() + t*b.real(),
        onet*a.imag() + t*b.imag()
    );
}


inline const complex& sum(const complex& c)
{
    return c;
}


template<class Cmpt> class Tensor;

//- No-op rotational transform for complex
inline complex transform(const Tensor<scalar>&, const complex c)
{
    return c;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline complex operator-(const complex& c)
{
    return complex(-c.re, -c.im);
}


inline complex operator+(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re + c2.re,
        c1.im + c2.im
    );
}


inline complex operator+(const complex& c, const scalar s)
{
    return complex(c.re + s, c.im);
}


inline complex operator+(const scalar s, const complex& c)
{
    return complex(c.re + s, c.im);
}


inline complex operator-(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re - c2.re,
        c1.im - c2.im
    );
}


inline complex operator-(const complex& c, const scalar s)
{
    return complex(c.re - s, c.im);
}


inline complex operator-(const scalar s, const complex& c)
{
    return complex(s - c.re, -c.im);
}


inline complex operator*(const complex& c1, const complex& c2)
{
    return complex
    (
        c1.re*c2.re - c1.im*c2.im,
        c1.im*c2.re + c1.re*c2.im
    );
}


inline complex operator*(const complex& c, const scalar s)
{
    return complex(s*c.re, s*c.im);
}


inline complex operator*(const scalar s, const complex& c)
{
    return complex(s*c.re, s*c.im);
}


inline complex operator/(const complex& c1, const complex& c2)
{
    const scalar sqrC2 = c2.magSqr();

    return complex
    (
        (c1.re*c2.re + c1.im*c2.im)/sqrC2,
        (c1.im*c2.re - c1.re*c2.im)/sqrC2
    );
}


inline complex operator/(const complex& c, const scalar s)
{
    return complex(c.re/s, c.im/s);
}


inline complex operator/(const scalar s, const complex& c)
{
    return complex(s)/c;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Complex transcendental functions
namespace Foam
{
    #define transFunc(func)                                                   \
    inline complex func(const complex& c)                                     \
    {                                                                         \
        return std:: func (static_cast<std::complex<scalar>>(c));             \
    }

transFunc(sqrt)
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

// Special treatment for pow()
inline complex pow(const complex& x, const complex& y)
{
    return std::pow
    (
        static_cast<std::complex<scalar>>(x),
        static_cast<std::complex<scalar>>(y)
    );
}


// Combinations of complex and integral/float
#define powFuncs(type2)                                                       \
    inline complex pow(const complex& x, const type2 y)                       \
    {                                                                         \
        return std::pow(static_cast<std::complex<scalar>>(x), y);             \
    }                                                                         \
                                                                              \
    inline complex pow(const type2 x, const complex& y)                       \
    {                                                                         \
        return std::pow                                                       \
        (                                                                     \
            static_cast<std::complex<scalar>>(x),                             \
            static_cast<std::complex<scalar>>(y)                              \
        );                                                                    \
    }

powFuncs(int)
powFuncs(long)
#if defined(__APPLE__) && WM_LABEL_SIZE == 64
powFuncs(int64_t)
#endif
powFuncs(float)
powFuncs(double)


inline complex pow3(const complex& c)
{
    return c*sqr(c);
}


inline complex pow4(const complex& c)
{
    return sqr(sqr(c));
}


inline complex pow5(const complex& c)
{
    return c*pow4(c);
}


inline complex pow6(const complex& c)
{
    return pow3(sqr(c));
}


inline complex pow025(const complex& c)
{
    return sqrt(sqrt(c));
}

} // End namespace Foam

#undef transFunc
#undef powFuncs

// ************************************************************************* //
