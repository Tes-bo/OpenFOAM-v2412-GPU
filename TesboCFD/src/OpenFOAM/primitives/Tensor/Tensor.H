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
    Foam::Tensor

Description
    A templated (3 x 3) tensor of objects of \<T\> derived from MatrixSpace.

SourceFiles
    TensorI.H

See also
    Foam::MatrixSpace
    Foam::Vector

\*---------------------------------------------------------------------------*/

#ifndef Foam_Tensor_H
#define Foam_Tensor_H

#include "MatrixSpace.H"
#include "Vector.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Cmpt> class SymmTensor;

/*---------------------------------------------------------------------------*\
                           Class Tensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class Tensor
:
    public MatrixSpace<Tensor<Cmpt>, Cmpt, 3, 3>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef Tensor<label> labelType;


    // Member Constants

        //- Rank of Tensor is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const Tensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YX, YY, YZ, ZX, ZY, ZZ };


    // Generated Methods

        //- Default construct
        Tensor() = default;

        //- Copy construct
        Tensor(const Tensor&) = default;

        //- Copy assignment
        Tensor& operator=(const Tensor&) = default;


    // Constructors

        //- Construct initialized to zero
        inline Tensor(const Foam::zero);

        //- Construct given MatrixSpace of the same rank
        template<class Cmpt2>
        inline Tensor(const MatrixSpace<Tensor<Cmpt2>, Cmpt2, 3, 3>& vs);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline Tensor(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>& vs);

        //- Construct given SphericalTensor
        inline Tensor(const SphericalTensor<Cmpt>& st);

        //- Construct given SymmTensor
        inline Tensor(const SymmTensor<Cmpt>& st);

        //- Construct given triad of row vectors,
        //- optionally treated as transposed (ie, column vectors)
        inline Tensor
        (
            const Vector<Vector<Cmpt>>& vecs,
            const bool transposed = false
        );

        //- Construct given the three row vectors,
        //- optionally treated as transposed (ie, column vectors)
        inline Tensor
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z,
            const bool transposed = false
        );

        //- Construct given the nine components
        inline Tensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
            const Cmpt tyx, const Cmpt tyy, const Cmpt tyz,
            const Cmpt tzx, const Cmpt tzy, const Cmpt tzz
        );

        //- Construct from a block of another matrix space
        template
        <
            template<class, direction, direction> class Block2,
            direction BRowStart,
            direction BColStart
        >
        Tensor
        (
            const Block2<Tensor<Cmpt>, BRowStart, BColStart>& block
        );

        //- Construct from Istream
        inline explicit Tensor(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& xx() const noexcept { return this->v_[XX]; }
        const Cmpt& xy() const noexcept { return this->v_[XY]; }
        const Cmpt& xz() const noexcept { return this->v_[XZ]; }
        const Cmpt& yx() const noexcept { return this->v_[YX]; }
        const Cmpt& yy() const noexcept { return this->v_[YY]; }
        const Cmpt& yz() const noexcept { return this->v_[YZ]; }
        const Cmpt& zx() const noexcept { return this->v_[ZX]; }
        const Cmpt& zy() const noexcept { return this->v_[ZY]; }
        const Cmpt& zz() const noexcept { return this->v_[ZZ]; }

        Cmpt& xx() noexcept { return this->v_[XX]; }
        Cmpt& xy() noexcept { return this->v_[XY]; }
        Cmpt& xz() noexcept { return this->v_[XZ]; }
        Cmpt& yx() noexcept { return this->v_[YX]; }
        Cmpt& yy() noexcept { return this->v_[YY]; }
        Cmpt& yz() noexcept { return this->v_[YZ]; }
        Cmpt& zx() noexcept { return this->v_[ZX]; }
        Cmpt& zy() noexcept { return this->v_[ZY]; }
        Cmpt& zz() noexcept { return this->v_[ZZ]; }


    // Column-vector access

        //- Extract vector for column 0
        inline Vector<Cmpt> cx() const;

        //- Extract vector for column 1
        inline Vector<Cmpt> cy() const;

        //- Extract vector for column 2
        inline Vector<Cmpt> cz() const;

        //- Extract vector for given column: compile-time check of index
        template<direction Idx>
        inline Vector<Cmpt> col() const;

        //- Extract vector for given column (0,1,2): runtime check of index
        inline Vector<Cmpt> col(const direction c) const;

        //- Set values of given column: compile-time check of index
        template<direction Idx>
        inline void col(const Vector<Cmpt>& v);

        //- Set values of given column (0,1,2): runtime check of index
        inline void col(const direction c, const Vector<Cmpt>& v);

        //- Set column values
        inline void cols
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z
        );


    // Row-vector access

        //- Extract vector for row 0
        inline Vector<Cmpt> x() const;

        //- Extract vector for row 1
        inline Vector<Cmpt> y() const;

        //- Extract vector for row 2
        inline Vector<Cmpt> z() const;

        //- Extract vector for given row: compile-time check of index
        template<direction Idx>
        inline Vector<Cmpt> row() const;

        //- Extract vector for given row (0,1,2): runtime check of index
        inline Vector<Cmpt> row(const direction r) const;

        //- Set values of given row: compile-time check of index
        template<direction Idx>
        inline void row(const Vector<Cmpt>& v);

        //- Set values of given row (0,1,2): runtime check of row
        inline void row(const direction r, const Vector<Cmpt>& v);

        //- Set row values
        inline void rows
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z
        );


    // Diagonal access and manipulation

        //- Extract the diagonal as a vector
        inline Vector<Cmpt> diag() const;

        //- Set values of the diagonal
        inline void diag(const Vector<Cmpt>& v);

        //- Add to the diagonal
        inline void addDiag(const Vector<Cmpt>& v);

        //- Subtract from the diagonal
        inline void subtractDiag(const Vector<Cmpt>& v);

        //- The L2-norm squared of the diagonal
        inline scalar diagSqr() const;


    // Characteristics

        //- Is identity tensor?
        inline bool is_identity(const scalar tol = ROOTVSMALL) const;


    // Tensor Operations

        //- Return non-Hermitian transpose
        inline Tensor<Cmpt> T() const;

        //- The determinate
        inline Cmpt det() const;

        //- The 2D determinant by excluding given direction
        inline Cmpt det2D(const direction excludeCmpt) const;

        //- Return adjunct matrix (transpose of cofactor matrix)
        inline Tensor<Cmpt> adjunct() const;

        //- Return cofactor matrix (transpose of adjunct matrix)
        inline Tensor<Cmpt> cof() const;

        //- Return 2D adjunct matrix by excluding given direction
        inline Tensor<Cmpt> adjunct2D(const direction excludeCmpt) const;

        //- Return inverse
        inline Tensor<Cmpt> inv() const;

        //- Return inverse, with (ad hoc) failsafe handling of 2D tensors
        inline Tensor<Cmpt> safeInv() const;

        //- Return inverse of 2D tensor (by excluding given direction)
        inline Tensor<Cmpt> inv2D(const direction excludeCmpt) const;

        //- Inner-product of this with another Tensor.
        inline Tensor<Cmpt> inner(const Tensor<Cmpt>& t2) const;

        //- Schur-product of this with another Tensor.
        inline Tensor<Cmpt> schur(const Tensor<Cmpt>& t2) const;


    // Member Operators

        //- Assign inner-product of this with another Tensor.
        inline void operator&=(const Tensor<Cmpt>& t);

        //- Inherit MatrixSpace assignment operators
        using Tensor::msType::operator=;

        //- Assign to an equivalent vector space
        template<class Cmpt2>
        inline void operator=(const VectorSpace<Tensor<Cmpt2>, Cmpt2, 9>&);

        //- Assign to a SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);

        //- Assign to a SymmTensor
        inline void operator=(const SymmTensor<Cmpt>&);

        //- Assign to a triad of row vectors
        inline void operator=(const Vector<Vector<Cmpt>>&);


    // Housekeeping

        //- Deprecated(2018-12) Return vector for given row (0,1)
        //  \deprecated(2018-12) use row() method
        Vector<Cmpt> vectorComponent(const direction cmpt) const
        {
            return row(cmpt);
        }
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<Tensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<Tensor<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<Tensor<Cmpt>> : is_contiguous_scalar<Cmpt> {};


template<class Cmpt>
class typeOfRank<Cmpt, 2>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfTranspose<Cmpt, Tensor<Cmpt>>
{
public:

    typedef Tensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSolve<Tensor<Cmpt>>
{
public:

    typedef Tensor<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Include inline implementations
#include "TensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
