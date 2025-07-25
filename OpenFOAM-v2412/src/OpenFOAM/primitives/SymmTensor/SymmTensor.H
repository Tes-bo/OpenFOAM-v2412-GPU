/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::SymmTensor

Description
    A templated (3 x 3) symmetric tensor of objects of \<T\>, effectively
    containing 6 elements, derived from VectorSpace.

SourceFiles
    SymmTensorI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_SymmTensor_H
#define Foam_SymmTensor_H

#include "Vector.H"
#include "SphericalTensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class SymmTensor Declaration
\*---------------------------------------------------------------------------*/

template<class Cmpt>
class SymmTensor
:
    public VectorSpace<SymmTensor<Cmpt>, Cmpt, 6>
{
public:

    // Typedefs

        //- Equivalent type of labels used for valid component indexing
        typedef SymmTensor<label> labelType;


    // Member Constants

        //- Rank of SymmTensor is 2
        static constexpr direction rank = 2;


    // Static Data Members

        static const SymmTensor I;


    //- Component labeling enumeration
    enum components { XX, XY, XZ, YY, YZ, ZZ };


    // Generated Methods

        //- Default construct
        SymmTensor() = default;

        //- Copy construct
        SymmTensor(const SymmTensor&) = default;

        //- Copy assignment
        SymmTensor& operator=(const SymmTensor&) = default;


    // Constructors

        //- Construct initialized to zero
        inline SymmTensor(const Foam::zero);

        //- Construct given VectorSpace of the same rank
        template<class Cmpt2>
        inline SymmTensor(const VectorSpace<SymmTensor<Cmpt2>, Cmpt2, 6>&);

        //- Construct given SphericalTensor
        inline SymmTensor(const SphericalTensor<Cmpt>&);

        //- Construct given the three row (or column) vectors
        inline SymmTensor
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z,
            const bool transposed = false  /* ignored */
        );

        //- Construct given the six components
        inline SymmTensor
        (
            const Cmpt txx, const Cmpt txy, const Cmpt txz,
                            const Cmpt tyy, const Cmpt tyz,
                                            const Cmpt tzz
        );

        //- Construct from Istream
        inline explicit SymmTensor(Istream& is);


    // Member Functions

    // Component Access

        const Cmpt& xx() const noexcept { return this->v_[XX]; }
        const Cmpt& xy() const noexcept { return this->v_[XY]; }
        const Cmpt& xz() const noexcept { return this->v_[XZ]; }
        const Cmpt& yx() const noexcept { return this->v_[XY]; }
        const Cmpt& yy() const noexcept { return this->v_[YY]; }
        const Cmpt& yz() const noexcept { return this->v_[YZ]; }
        const Cmpt& zx() const noexcept { return this->v_[XZ]; }
        const Cmpt& zy() const noexcept { return this->v_[YZ]; }
        const Cmpt& zz() const noexcept { return this->v_[ZZ]; }

        Cmpt& xx() noexcept { return this->v_[XX]; }
        Cmpt& xy() noexcept { return this->v_[XY]; }
        Cmpt& xz() noexcept { return this->v_[XZ]; }
        Cmpt& yx() noexcept { return this->v_[XY]; }
        Cmpt& yy() noexcept { return this->v_[YY]; }
        Cmpt& yz() noexcept { return this->v_[YZ]; }
        Cmpt& zx() noexcept { return this->v_[XZ]; }
        Cmpt& zy() noexcept { return this->v_[YZ]; }
        Cmpt& zz() noexcept { return this->v_[ZZ]; }


    // Column-vector access

        //- Extract vector for column 0
        Vector<Cmpt> cx() const { return this->x(); }

        //- Extract vector for column 1
        Vector<Cmpt> cy() const { return this->y(); }

        //- Extract vector for column 2
        Vector<Cmpt> cz() const { return this->z(); }

        //- Extract vector for given column: compile-time check of index
        template<direction Idx>
        Vector<Cmpt> col() const { return this->row<Idx>(); }

        //- Extract vector for given column (0,1,2): runtime check of index
        Vector<Cmpt> col(const direction c) const { return this->row(c); }

        //- Set values of given column (0,1,2): runtime check of index
        void col(const direction c, const Vector<Cmpt>& v) { this->row(c, v); }

        //- Set column values
        void cols
        (
            const Vector<Cmpt>& x,
            const Vector<Cmpt>& y,
            const Vector<Cmpt>& z
        )
        {
            this->rows(x, y, z);
        }


    // Row-vector access

        //- Extract vector for row 0
        inline Vector<Cmpt> x() const;

        //- Extract vector for row 1
        inline Vector<Cmpt> y() const;

        //- Extract vector for row 2
        inline Vector<Cmpt> z() const;

        //- Extract vector for given row: compile-time check of index
        template<direction Row>
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


    // Tensor Operations

        //- Return non-Hermitian transpose
        const SymmTensor<Cmpt>& T() const noexcept { return *this; }

        //- The determinate
        inline Cmpt det() const;

        //- The 2D determinant by excluding given direction
        inline Cmpt det2D(const direction excludeCmpt) const;

        //- Return adjunct matrix (transpose of cofactor matrix)
        inline SymmTensor<Cmpt> adjunct() const;

        //- Return cofactor matrix (transpose of adjunct matrix)
        inline SymmTensor<Cmpt> cof() const;

        //- Return 2D adjunct matrix by excluding given direction
        inline SymmTensor<Cmpt> adjunct2D(const direction excludeCmpt) const;

        //- Return inverse
        inline SymmTensor<Cmpt> inv() const;

        //- Return inverse, with (ad hoc) failsafe handling of 2D tensors
        inline SymmTensor<Cmpt> safeInv() const;

        //- Return inverse of 2D tensor (by excluding given direction)
        inline SymmTensor<Cmpt> inv2D(const direction excludeCmpt) const;


    // Member Operators

        //- Inherit VectorSpace assignment operators
        using SymmTensor::vsType::operator=;

        //- Assign to given SphericalTensor
        inline void operator=(const SphericalTensor<Cmpt>&);
};


// * * * * * * * * * * * * * * * * * Traits  * * * * * * * * * * * * * * * * //

//- Data are contiguous if component type is contiguous
template<class Cmpt>
struct is_contiguous<SymmTensor<Cmpt>> : is_contiguous<Cmpt> {};

//- Data are contiguous label if component type is label
template<class Cmpt>
struct is_contiguous_label<SymmTensor<Cmpt>> : is_contiguous_label<Cmpt> {};

//- Data are contiguous scalar if component type is scalar
template<class Cmpt>
struct is_contiguous_scalar<SymmTensor<Cmpt>> : is_contiguous_scalar<Cmpt> {};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 2>
{
public:

    typedef SymmTensor<Cmpt> type;
};


template<class Cmpt>
class typeOfSolve<SymmTensor<Cmpt>>
{
public:

    typedef SymmTensor<solveScalar> type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SymmTensorI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
