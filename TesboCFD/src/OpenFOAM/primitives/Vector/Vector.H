/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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
    Foam::Vector

Description
    Templated 3D Vector derived from VectorSpace adding construction from
    3 components, element access using x(), y() and z() member functions and
    the inner-product (dot-product) and cross-product operators.

    A centre() member function which returns the Vector for which it is called
    is defined so that point which is a typedef to Vector\<scalar\> behaves as
    other shapes in the shape hierarchy.

SourceFiles
    VectorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Vector_H
#define Foam_Vector_H

#include "contiguous.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class T> class List;

/*---------------------------------------------------------------------------*\
                           Class Vector Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Vector
:
    public VectorSpace<Vector<Cmpt>, Cmpt, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Vector<label> labelType;


    // Member Constants

        //- Rank of Vector is 1
        static constexpr direction rank = 1;


    //- Component labeling enumeration
    enum components { X, Y, Z };


    // Generated Methods

        //- Default construct
        Vector() = default;

        //- Copy construct
        Vector(const Vector&) = default;

        //- Copy assignment
        Vector& operator=(const Vector&) = default;


    // Constructors

        //- Construct initialized to zero
        inline Vector(const Foam::zero);

        //- Copy construct from VectorSpace of the same rank
        template<class Cmpt2>
        inline Vector(const VectorSpace<Vector<Cmpt2>, Cmpt2, 3>& vs);

        //- Construct from three components
        inline Vector(const Cmpt& vx, const Cmpt& vy, const Cmpt& vz);

        //- Construct from Istream
        inline explicit Vector(Istream& is);


    // Member Functions

    // Component Access

        //- Access to the vector x component
        const Cmpt& x() const noexcept { return this->v_[X]; }

        //- Access to the vector y component
        const Cmpt& y() const noexcept { return this->v_[Y]; }

        //- Access to the vector z component
        const Cmpt& z() const noexcept { return this->v_[Z]; }

        //- Access to the vector x component
        Cmpt& x() noexcept { return this->v_[X]; }

        //- Access to the vector y component
        Cmpt& y() noexcept { return this->v_[Y]; }

        //- Access to the vector z component
        Cmpt& z() noexcept { return this->v_[Z]; }


    // Vector Operations

        //- Return \c this (for point which is a typedef to Vector\<scalar\>)
        inline const Vector<Cmpt>& centre
        (
            const Foam::UList<Vector<Cmpt>>&  /* (unused) */
        ) const noexcept;

        //- The length (L2-norm) of the vector
        inline scalar mag() const;

        //- The length (L2-norm) squared of the vector.
        inline scalar magSqr() const;

        //- The L2-norm distance from another vector.
        //- The mag() of the difference.
        inline scalar dist(const Vector<Cmpt>& v2) const;

        //- The L2-norm distance squared from another vector.
        //- The magSqr() of the difference.
        inline scalar distSqr(const Vector<Cmpt>& v2) const;

        //- Inplace normalise the vector by its magnitude
        //  For small magnitudes (less than ROOTVSMALL) set to zero.
        //  Will not be particularly useful for a vector of labels
        inline Vector<Cmpt>& normalise(const scalar tol = ROOTVSMALL);

        //- Inplace removal of components that are collinear to the given
        //- unit vector.
        inline Vector<Cmpt>& removeCollinear(const Vector<Cmpt>& unitVec);

        //- Scalar-product of \c this with another Vector.
        inline Cmpt inner(const Vector<Cmpt>& v2) const;

        //- Cross-product of \c this with another Vector.
        inline Vector<Cmpt> cross(const Vector<Cmpt>& v2) const;


    // Comparison Operations

        //- Lexicographically compare \em a and \em b with order (x:y:z)
        static inline bool less_xyz
        (
            const Vector<Cmpt>& a,
            const Vector<Cmpt>& b
        );

        //- Lexicographically compare \em a and \em b with order (y:z:x)
        static inline bool less_yzx
        (
            const Vector<Cmpt>& a,
            const Vector<Cmpt>& b
        );

        //- Lexicographically compare \em a and \em b with order (z:x:y)
        static inline bool less_zxy
        (
            const Vector<Cmpt>& a,
            const Vector<Cmpt>& b
        );
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Vector<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Vector<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Vector<Cmpt>> : is_contiguous_scalar<Cmpt> {};


template<class Cmpt>
class typeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 1>
{
public:

    typedef Vector<Cmpt> type;
};


template<class Cmpt>
class typeOfSolve<Vector<Cmpt>>
{
public:

    typedef Vector<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "VectorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
