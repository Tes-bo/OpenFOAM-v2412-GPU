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

Class
    Foam::complex

Description
    A complex number, similar to the C++ complex type.

SourceFiles
    complexI.H
    complex.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_primitives_complex_H
#define Foam_primitives_complex_H

#include "scalar.H"
#include "word.H"
#include "zero.H"
#include "contiguous.H"
#include <complex>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class complex;

inline complex operator-(const complex&);
inline complex operator+(const complex&, const complex&);
inline complex operator+(const complex&, const scalar);
inline complex operator+(const scalar, const complex&);
inline complex operator-(const complex&, const complex&);
inline complex operator-(const complex&, const scalar);
inline complex operator-(const scalar, const complex&);
inline complex operator*(const complex&, const complex&);
inline complex operator*(const complex&, const scalar);
inline complex operator*(const scalar, const complex&);
inline complex operator/(const complex&, const complex&);
inline complex operator/(const complex&, const scalar);
inline complex operator/(const scalar, const complex&);


/*---------------------------------------------------------------------------*\
                           Class complex Declaration
\*---------------------------------------------------------------------------*/

class complex
{
    // Private Data

        //- Real and imaginary parts
        scalar re, im;


public:

    // Generated Methods

        //- Copy construct
        complex(const complex&) noexcept = default;

        //- Copy assignment
        complex& operator=(const complex&) noexcept = default;

        //- Move construct
        complex(complex&&) noexcept = default;

        //- Move assignment
        complex& operator=(complex&&) noexcept = default;


    // Constructors

        //- Default construct, as zero-initialized
        inline constexpr complex() noexcept;

        //- Construct zero-initialized from zero class
        inline constexpr complex(const Foam::zero) noexcept;

        //- Construct from real component
        inline explicit constexpr complex(const scalar r) noexcept;

        //- Construct from real and imaginary parts
        inline constexpr complex(const scalar r, const scalar i) noexcept;

        //- Implicit construct from std::complex
        inline complex(const std::complex<float>& c);

        //- Implicit construct from std::complex
        inline complex(const std::complex<double>& c);

        //- Construct from Istream
        explicit complex(Istream& is);


    // Member Functions

    // STL getter/setter

        //- Real part of complex number - STL naming
        constexpr scalar real() const noexcept { return re; }

        //- Imaginary part of complex number - STL naming
        constexpr scalar imag() const noexcept { return im; }

        //- Set real part of complex number - STL naming
        void real(scalar val) noexcept { re = val; }

        //- Set imaginary part of complex number - STL naming
        void imag(scalar val) noexcept { im = val; }


    // Methods

        //- Complex conjugate
        inline complex conjugate() const;

        //- The magnitude (L2-norm) of complex.
        //- Called magnitude() instead mag(), which looks too much like imag()
        inline scalar magnitude() const;

        //- The L2-norm squared of complex
        inline scalar magSqr() const;

        //- The sum of real/imag components
        inline scalar cmptSum() const noexcept;


    // Member Operators

        //- Implicit conversion to std::complex
        operator std::complex<scalar>() const
        {
            return std::complex<scalar>(re, im);
        }


        //- Assign zero
        inline void operator=(const Foam::zero);

        //- Assign scalar (imag = zero)
        inline void operator=(const scalar s);

        inline void operator+=(const complex& c);
        inline void operator+=(const scalar s);

        inline void operator-=(const complex& c);
        inline void operator-=(const scalar s);

        inline void operator*=(const complex& c);
        inline void operator*=(const scalar s);

        inline void operator/=(const complex& c);
        inline void operator/=(const scalar s);

        inline bool operator==(const complex& c) const;
        inline bool operator!=(const complex& c) const;


    // Friend Operators

        friend complex operator-(const complex& c);

        friend complex operator+(const complex& c1, const complex& c2);
        friend complex operator+(const complex& c, const scalar s);
        friend complex operator+(const scalar s, const complex& c);

        friend complex operator-(const complex& c1, const complex& c2);
        friend complex operator-(const complex& c, const scalar s);
        friend complex operator-(const scalar s, const complex& c);

        friend complex operator*(const complex& c1, const complex& c2);
        friend complex operator*(const complex& c, const scalar s);
        friend complex operator*(const scalar s, const complex& c);

        friend complex operator/(const complex& c1, const complex& c2);
        friend complex operator/(const complex& c, const scalar s);
        friend complex operator/(const scalar s, const complex& c);


    // Housekeeping

        //- Get real part of complex number. Same as real()
        scalar Re() const noexcept { return re; }

        //- Get imaginary part of complex number. Same as imag()
        scalar Im() const noexcept { return im; }

        //- Non-const access to real part. Prefer real() setter method
        scalar& Re() noexcept { return re; }

        //- Non-const access to imaginary part. Prefer imag() setter method
        scalar& Im() noexcept { return im; }
};


/*---------------------------------------------------------------------------*\
                        Specialization pTraits<complex>
\*---------------------------------------------------------------------------*/

//- The underlying component data type for complex is scalar.
//  The regular pTraits<T>:cmptType as complex is currently (2023-11)
//  likely not quite correct (issue #3018)
template<>
struct pTraits_cmptType<complex> { typedef scalar type; };

//- A complex has two scalar components
template<>
struct pTraits_nComponents<complex>
:
    std::integral_constant<Foam::direction, 2>
{};


//- Template specialisation for pTraits<complex>
template<>
class pTraits<complex>
{
    complex p_;

public:

    // Typedefs

        //- Component type
        typedef complex cmptType;

        //- Magnitude type
        typedef scalar magType;

        //- Equivalent type of labels used for valid component indexing
        typedef label labelType;


    // Member Constants

        //- Dimensionality of space
        static constexpr direction dim = 3;

        //- Rank of complex is 0
        static constexpr direction rank = 0;

        //- Number of components in complex is 2
        static constexpr direction nComponents = 2;


    // Static Data Members

        static const char* const typeName;
        static const char* const componentNames[];

        static const complex zero;      //!< complex (0,0)
        static const complex one;       //!< complex (1,0)
        static const complex min;       //!< complex (-VGREAT,-VGREAT)
        static const complex max;       //!< complex (VGREAT,VGREAT)
        static const complex rootMin;   //!< complex (-ROOTVGREAT, -ROOTVGREAT)
        static const complex rootMax;   //!< complex (ROOTVGREAT, ROOTVGREAT)


    // Constructors

        //- Copy construct from primitive
        explicit pTraits(const complex& val) noexcept
        :
            p_(val)
        {}


        //- Read construct from Istream
        explicit pTraits(Istream& is);


    // Member Functions

        //- Return the value
        operator complex() const noexcept { return p_; }

        //- Access the value
        operator complex&() noexcept { return p_; }
};


/*---------------------------------------------------------------------------*\
                               Namespace Detail
\*---------------------------------------------------------------------------*/

namespace Detail
{
    // Helper functions for complex, in Detail namespace to avoid possible
    // name collisions (could change in the future)

    //- The 'conjugate' of non-complex returns itself (pass-through)
    //- it does not return a complex!
    template<class T>
    typename std::enable_if
    <
        !std::is_same<complex, T>::value,
        const T&
    >::type conj(const T& val)
    {
        return val;
    }

    //- The conjugate of a complex number
    template<class T>
    typename std::enable_if
    <
        std::is_same<complex, T>::value,
        complex
    >::type conj(const T& val)
    {
        return val.conjugate();
    }

} // End namespace Detail


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Contiguous data for complex
template<> struct is_contiguous<complex> : std::true_type {};

//- Contiguous scalar data for complex
template<> struct is_contiguous_scalar<complex> : std::true_type {};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

Istream& operator>>(Istream& is, complex& c);
Ostream& operator<<(Ostream& os, const complex& c);

//- Complex conjugate
inline complex operator~(const complex& c);


// * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * * //

//- Return string representation of complex
word name(const complex& c);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "complexI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
