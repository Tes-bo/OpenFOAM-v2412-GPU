/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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
    Foam::Tensor2D

Description
    A templated (2 x 2) tensor of objects of \<T\> derived from VectorSpace.

SourceFiles
    Tensor2DI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_Tensor2D_H
#define Foam_Tensor2D_H

#include "Vector2D.H"
#include "SphericalTensor2D.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Cmpt> class SymmTensor2D;

/*---------------------------------------------------------------------------*\
                           Class Tensor2D Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor2D
:
    public VectorSpace<Tensor2D<Cmpt>, Cmpt, 4>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Tensor2D<label> labelType;


    // Member Constants

        //- Rank of Tensor2D is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const Tensor2D I;


    //- Component labeling enumeration
    enum components { XX, XY, YX, YY };


    // Generated Methods

        //- Default construct
        Tensor2D() = default;

        //- Copy construct
        Tensor2D(const Tensor2D&) = default;

        //- Copy assignment
        Tensor2D& operator=(const Tensor2D&) = default;


    // Constructors

        //- Construct initialized to zero
        inline Tensor2D(const Foam::zero);

        //- Construct given VectorSpace
        inline Tensor2D(const VectorSpace<Tensor2D<Cmpt>, Cmpt, 4>& vs);

        //- Construct given SymmTensor2D
        inline Tensor2D(const SymmTensor2D<Cmpt>& st);

        //- Construct given SphericalTensor2D
        inline Tensor2D(const SphericalTensor2D<Cmpt>& st);

        //- Construct given the two row vectors
        inline Tensor2D
        (
            const Vector2D<Cmpt>& x,
            const Vector2D<Cmpt>& y
        );

        //- Construct given the four components
        inline Tensor2D
        (
            const Cmpt txx, const Cmpt txy,
            const Cmpt tyx, const Cmpt tyy
        );

        //- Construct from Istream
        inline explicit Tensor2D(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& xx() const noexcept { return this->v_[XX]; }
        const Cmpt& xy() const noexcept { return this->v_[XY]; }
        const Cmpt& yx() const noexcept { return this->v_[YX]; }
        const Cmpt& yy() const noexcept { return this->v_[YY]; }

        Cmpt& xx() noexcept { return this->v_[XX]; }
        Cmpt& xy() noexcept { return this->v_[XY]; }
        Cmpt& yx() noexcept { return this->v_[YX]; }
        Cmpt& yy() noexcept { return this->v_[YY]; }


    // Column-vector access

        //- Extract vector for column 0
        inline Vector2D<Cmpt> cx() const;

        //- Extract vector for column 1
        inline Vector2D<Cmpt> cy() const;

        //- Extract vector for given column: compile-time check of index
        template<direction Idx>
        inline Vector2D<Cmpt> col() const;

        //- Extract vector for given column (0,1): runtime check of index
        inline Vector2D<Cmpt> col(const direction c) const;

        //- Set values of given column: compile-time check of index
        template<direction Idx>
        inline void col(const Vector2D<Cmpt>& v);

        //- Set values of given column (0,1): runtime check of index
        inline void col(const direction c, const Vector2D<Cmpt>& v);

        //- Set column values
        inline void cols(const Vector2D<Cmpt>& x, const Vector2D<Cmpt>& y);


    // Row-vector access

        //- Extract vector for row 0
        inline Vector2D<Cmpt> x() const;

        //- Extract vector for row 1
        inline Vector2D<Cmpt> y() const;

        //- Extract vector for given row: compile-time check of index
        template<direction Idx>
        inline Vector2D<Cmpt> row() const;

        //- Extract vector for given row (0,1): runtime check of index
        inline Vector2D<Cmpt> row(const direction r) const;

        //- Set values of given row: compile-time check of index
        template<direction Idx>
        inline void row(const Vector2D<Cmpt>& v);

        //- Set values of given row (0,1): compile-time check of index
        inline void row(const direction r, const Vector2D<Cmpt>& v);

        //- Set row values
        inline void rows(const Vector2D<Cmpt>& x, const Vector2D<Cmpt>& y);


    // Diagonal access and manipulation

        //- Extract the diagonal as a vector
        inline Vector2D<Cmpt> diag() const;

        //- Set values of the diagonal
        inline void diag(const Vector2D<Cmpt>& v);

        //- The L2-norm squared of the diagonal
        inline scalar diagSqr() const;


    // Tensor Operations

        //- Return non-Hermitian transpose
        inline Tensor2D<Cmpt> T() const;

        //- The determinate
        inline Cmpt det() const;

        //- Return adjunct matrix (transpose of cofactor matrix)
        inline Tensor2D<Cmpt> adjunct() const;

        //- Return cofactor matrix (transpose of adjunct matrix)
        inline Tensor2D<Cmpt> cof() const;

        //- Return inverse
        inline Tensor2D<Cmpt> inv() const;

        //- Inner-product of this with another Tensor2D.
        inline Tensor2D<Cmpt> inner(const Tensor2D<Cmpt>& t2) const;

        //- Schur-product of this with another Tensor2D.
        inline Tensor2D<Cmpt> schur(const Tensor2D<Cmpt>& t2) const;


    // Member Operators

        //- Copy assign from SymmTensor2D
        inline void operator=(const SymmTensor2D<Cmpt>&);

        //- Copy assign from SphericalTensor2D
        inline void operator=(const SphericalTensor2D<Cmpt>&);


    // Housekeeping

        //- Deprecated(2018-12) Return vector for given row (0,1)
        //  \deprecated(2018-12) use row() method
        Vector2D<Cmpt> vectorComponent(const direction cmpt) const
        {
            return row(cmpt);
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Tensor2D<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Tensor2D<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Tensor2D<Cmpt>> : is_contiguous_scalar<Cmpt> {};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "Tensor2DI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
