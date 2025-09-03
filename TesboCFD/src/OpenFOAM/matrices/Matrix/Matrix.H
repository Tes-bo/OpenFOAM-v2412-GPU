/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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
    Foam::Matrix

Description
    A templated (m x n) matrix of objects of \<T\>.
    The layout is (mRows x nCols) - row-major order:

    \verbatim
        (0,0)
          +---> j [nCols]
          |
          |
          v
          i [mRows]
    \endverbatim

SourceFiles
    Matrix.C
    MatrixI.H
    MatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_Matrix_H
#define Foam_Matrix_H

#include "uLabel.H"
#include "Pair.H"
#include "Field.H"
#include "autoPtr.H"
#include "complex.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
template<class Form, class Type> class Matrix;
template<class MatrixType> class ConstMatrixBlock;
template<class MatrixType> class MatrixBlock;

template<class Form, class Type>
Istream& operator>>(Istream& is, Matrix<Form, Type>& mat);

template<class Form, class Type>
Ostream& operator<<(Ostream& os, const Matrix<Form, Type>& mat);


/*---------------------------------------------------------------------------*\
                           Class Matrix Declaration
\*---------------------------------------------------------------------------*/

template<class Form, class Type>
class Matrix
{
    // Private Data

        //- Number of rows and columns in Matrix
        label mRows_, nCols_;

        //- Vector of values of type Type
        Type* __restrict__ v_;


    // Private Member Functions

        //- Allocate storage for the contents
        inline void doAlloc();

        //- Right-multiply Matrix by a column vector (A * x)
        template<class ListType>
        tmp<Field<Type>> AmulImpl(const ListType& x) const;

        //- Left-multiply Matrix by a row vector (x * A)
        template<class ListType>
        tmp<Field<Type>> TmulImpl(const ListType& x) const;


public:

    //- Matrix type
    typedef Matrix<Form, Type> mType;

    //- Component type
    typedef Type cmptType;


    // Static Member Functions

        //- Return a null Matrix (reference to a nullObject).
        //- Behaves like a empty Matrix.
        static const Matrix<Form, Type>& null() noexcept
        {
            return NullObjectRef<Matrix<Form, Type>>();
        }


    // Iterators

        //- Random access iterator for traversing a Matrix
        typedef Type* iterator;

        //- Random access iterator for traversing a Matrix
        typedef const Type* const_iterator;


    // Constructors

        //- Default construct (empty matrix)
        inline Matrix() noexcept;

        //- Construct given number of rows/columns, uninitialised content
        Matrix(const label m, const label n);

        //- Construct with given number of rows/columns
        //- initializing all elements to zero
        Matrix(const label m, const label n, const Foam::zero);

        //- Construct with given number of rows/columns
        //- initializing all elements to the given value
        Matrix(const label m, const label n, const Type& val);

        //- Construct given number of rows/columns
        inline explicit Matrix(const labelPair& dims);

        //- Construct given number of rows/columns
        //- initializing all elements to zero
        inline Matrix(const labelPair& dims, const Foam::zero);

        //- Construct with given number of rows/columns
        //- initializing all elements to the given value
        inline Matrix(const labelPair& dims, const Type& val);

        //- Copy construct
        Matrix(const Matrix<Form, Type>& mat);

        //- Move construct
        Matrix(Matrix<Form, Type>&& mat);

        //- Copy constructor from Matrix of a different form
        template<class Form2>
        explicit Matrix(const Matrix<Form2, Type>& mat);

        //- Construct from a block of another Matrix
        template<class MatrixType>
        Matrix(const ConstMatrixBlock<MatrixType>& Mb);

        //- Construct from a block of another Matrix
        template<class MatrixType>
        Matrix(const MatrixBlock<MatrixType>& Mb);

        //- Construct from Istream.
        explicit Matrix(Istream& is);

        //- Clone
        inline autoPtr<mType> clone() const;


    //- Destructor
    ~Matrix();


    // Member Functions

    // Access

        //- The number of rows
        label mRows() const noexcept { return mRows_; }

        //- The number of rows
        label nRows() const noexcept { return mRows_; }

        //- The number of rows
        label m() const noexcept { return mRows_; }

        //- The number of columns
        label nCols() const noexcept { return nCols_; }

        //- The number of columns
        label n() const noexcept { return nCols_; }

        //- Return true if Matrix is empty (i.e., size() is zero)
        inline bool empty() const noexcept;

        //- The number of elements in Matrix (m*n)
        inline label size() const;

        //- Return row/column sizes
        inline labelPair sizes() const;

        //- Return const pointer to the first data element, which can also
        //- be used to address into Matrix contents
        inline const Type* cdata() const noexcept;

        //- Return pointer to the first data element, which can also
        //- be used to address into Matrix contents
        inline Type* data() noexcept;

        //- Return pointer to the underlying array serving as data storage,
        //- reinterpreted as byte data.
        //  \note Only meaningful for contiguous data
        inline const char* cdata_bytes() const noexcept;

        //- Return pointer to the underlying array serving as data storage,
        //- reinterpreted as byte data.
        //  \note Only meaningful for contiguous data
        inline char* data_bytes() noexcept;

        //- Number of contiguous bytes for the Matrix data,
        //- no runtime check that the type is actually contiguous
        //  \note Only meaningful for contiguous data
        inline std::streamsize size_bytes() const noexcept;

        //- Number of contiguous bytes for the Matrix data,
        //- runtime FatalError if type is not contiguous
        std::streamsize byteSize() const;

        //- Return const pointer to data in the specified row
        //  Subscript checking only with FULLDEBUG
        inline const Type* rowData(const label irow) const;

        //- Return pointer to data in the specified row
        //  Subscript checking only with FULLDEBUG
        inline Type* rowData(const label irow);

        //- Linear addressing const element access
        //  Subscript checking only with FULLDEBUG
        inline const Type& at(const label idx) const;

        //- Linear addressing element access
        //  Subscript checking only with FULLDEBUG
        inline Type& at(const label idx);

    // Block Access (const)

        //- Return const column or column's subset of Matrix
        //  Return entire column by its index: M.subColumn(a);
        //  Return subset of a column starting from rowIndex: M.subColumn(a, b);
        //  Return subset of a column starting from rowIndex with szRows elems:
        //      M.subColumn(a, b, c);
        inline ConstMatrixBlock<mType> subColumn
        (
            const label colIndex,
            const label rowIndex = 0,
            label len = -1
        ) const;

        //- Return const row or const row's subset of Matrix
        //  Return entire row by its index: M.subRow(a);
        //  Return subset of a row starting from columnIndex: M.subRow(a,b);
        //  Return subset of a row starting from columnIndex with szCols elems:
        //      M.subRow(a, b, c);
        inline ConstMatrixBlock<mType> subRow
        (
            const label rowIndex,
            const label colIndex = 0,
            label len = -1
        ) const;

        //- Return const sub-block of Matrix
        //  Sub-block starting at columnIndex & rowIndex indices
        inline ConstMatrixBlock<mType> subMatrix
        (
            const label rowIndex,
            const label colIndex,
            label szRows = -1,
            label szCols = -1
        ) const;

        //- Access Field as a ConstMatrixBlock
        template<class VectorSpace>
        inline ConstMatrixBlock<mType> block
        (
            const label rowIndex,
            const label colIndex
        ) const;

    // Block Access (non-const)

        //- Return column or column's subset of Matrix
        inline MatrixBlock<mType> subColumn
        (
            const label colIndex,
            const label rowIndex = 0,
            label len = -1
        );

        //- Return row or row's subset of Matrix
        inline MatrixBlock<mType> subRow
        (
            const label rowIndex,
            const label colIndex = 0,
            label len = -1
        );

        //- Return sub-block of Matrix
        inline MatrixBlock<mType> subMatrix
        (
            const label rowIndex,
            const label colIndex,
            label szRows = -1,
            label szCols = -1
        );

        //- Access Field as a MatrixBlock
        template<class VectorSpace>
        inline MatrixBlock<mType> block
        (
            const label rowIndex,
            const label colIndex
        );

    // Check

        //- Check index i is within valid range [0, m)
        inline void checki(const label irow) const;

        //- Check index j is within valid range [0, n)
        inline void checkj(const label jcol) const;

        //- Check that dimensions are positive, non-zero
        inline void checkSize() const;

        //- True if all entries have identical values, and Matrix is non-empty
        inline bool uniform() const;

    // Edit

        //- Clear Matrix, i.e. set sizes to zero
        void clear();

        //- Release storage management of Matrix contents by transferring
        //- management to a List
        List<Type> release();

        //- Swap contents
        void swap(Matrix<Form, Type>& mat);

        //- Transfer the contents of the argument Matrix into this Matrix
        //- and annul the argument Matrix
        void transfer(Matrix<Form, Type>& mat);

        //- Change Matrix dimensions, preserving the elements
        void resize(const label m, const label n);

        //- Change Matrix dimensions \em without preserving existing content
        void resize_nocopy(const label mrow, const label ncol);

        //- Change Matrix dimensions, preserving the elements
        inline void setSize(const label m, const label n);

        //- Resize Matrix without reallocating storage (unsafe)
        inline void shallowResize(const label m, const label n);

        //- Round elements with magnitude smaller than tol (SMALL) to zero
        void round(const scalar tol = SMALL);


    // Operations

        //- Return conjugate transpose of Matrix
        Form T() const;

        //- Return non-conjugate transpose of Matrix
        Form transpose() const;

        //- Right-multiply Matrix by a column vector (A * x)
        inline tmp<Field<Type>> Amul(const UList<Type>& x) const;

        //- Right-multiply Matrix by a column vector (A * x)
        template<class Addr>
        inline tmp<Field<Type>> Amul
        (
            const IndirectListBase<Type, Addr>& x
        ) const;

        //- Left-multiply Matrix by a row vector (x * A)
        inline tmp<Field<Type>> Tmul(const UList<Type>& x) const;

        //- Left-multiply Matrix by a row vector (x * A)
        template<class Addr>
        inline tmp<Field<Type>> Tmul
        (
            const IndirectListBase<Type, Addr>& x
        ) const;

        //- Extract the diagonal elements. Method may change in the future.
        List<Type> diag() const;

        //- Assign diagonal of Matrix
        void diag(const UList<Type>& list);

        //- Return the trace
        Type trace() const;

        //- Return L2-Norm of chosen column
        //  Optional without sqrt for parallel usage.
        scalar columnNorm(const label colIndex, const bool noSqrt=false) const;

        //- Return Frobenius norm of Matrix
        //  Optional without sqrt for parallel usage.
        scalar norm(const bool noSqrt=false) const;


    // Member Operators

        //- Return const pointer to data in the specified row - rowData().
        //  Subscript checking only with FULLDEBUG
        inline const Type* operator[](const label irow) const;

        //- Return pointer to data in the specified row - rowData().
        //  Subscript checking only with FULLDEBUG
        inline Type* operator[](const label irow);

        //- (i, j) const element access operator
        //  Subscript checking only with FULLDEBUG
        inline const Type& operator()(const label irow, const label jcol) const;

        //- (i, j) element access operator
        //  Subscript checking only with FULLDEBUG
        inline Type& operator()(const label irow, const label jcol);

        //- Copy assignment. Takes linear time.
        void operator=(const Matrix<Form, Type>& mat);

        //- Move assignment
        void operator=(Matrix<Form, Type>&& mat);

        //- Assignment to a block of another Matrix
        template<class MatrixType>
        void operator=(const ConstMatrixBlock<MatrixType>& Mb);

        //- Assignment to a block of another Matrix
        template<class MatrixType>
        void operator=(const MatrixBlock<MatrixType>& Mb);

        //- Assignment of all elements to zero
        void operator=(const Foam::zero);

        //- Assignment of all elements to the given value
        void operator=(const Type& val);

        //- Matrix addition
        void operator+=(const Matrix<Form, Type>& other);

        //- Matrix subtraction
        void operator-=(const Matrix<Form, Type>& other);

        //- Matrix scalar addition
        void operator+=(const Type& s);

        //- Matrix scalar subtraction
        void operator-=(const Type& s);

        //- Matrix scalar multiplication
        void operator*=(const Type& s);

        //- Matrix scalar division
        void operator/=(const Type& s);


    // Random Access Iterator (non-const)

        //- Return an iterator to begin traversing a Matrix
        inline iterator begin() noexcept;

        //- Return an iterator to end traversing a Matrix
        inline iterator end() noexcept;


    // Random Access Iterator (const)

        //- Return const_iterator to begin traversing a constant Matrix
        inline const_iterator cbegin() const noexcept;

        //- Return const_iterator to end traversing a constant Matrix
        inline const_iterator cend() const noexcept;

        //- Return const_iterator to begin traversing a constant Matrix
        inline const_iterator begin() const noexcept;

        //- Return const_iterator to end traversing a constant Matrix
        inline const_iterator end() const noexcept;


    // IO

        //- Read Matrix from Istream, discarding existing contents.
        bool readMatrix(Istream& is);

        //- Write Matrix, with line-breaks in ASCII when length exceeds
        //- shortLen.
        //  Using '0' suppresses line-breaks entirely.
        Ostream& writeMatrix(Ostream& os, const label shortLen=0) const;


    // Housekeeping

        //- Deprecated(2019-04) raw data pointer, const access
        //  \deprecated(2019-04) - use cdata() method
        FOAM_DEPRECATED_FOR(2019-04, "cdata() method")
        const Type* v() const
        {
            return v_;
        }

        //- Deprecated(2019-04) raw data pointer, non-const access
        //  \deprecated(2019-04) - use data() method
        FOAM_DEPRECATED_FOR(2019-04, "data() method")
        Type* v()
        {
            return v_;
        }

        //- Deprecated(2019-04) - use subMatrix()
        //  \deprecated(2019-04) - use subMatrix()
        FOAM_DEPRECATED_FOR(2019-04, "subMatrix() method")
        ConstMatrixBlock<mType> block
        (
            const label m,
            const label n,
            const label mStart,
            const label nStart
        ) const
        {
            return this->subMatrix(mStart, nStart, m, n);
        }

        //- Deprecated(2019-04) - use subMatrix()
        //  \deprecated(2019-04) - use subMatrix()
        FOAM_DEPRECATED_FOR(2019-04, "subMatrix() method")
        MatrixBlock<mType> block
        (
            const label m,
            const label n,
            const label mStart,
            const label nStart
        )
        {
            return this->subMatrix(mStart, nStart, m, n);
        }


        //- Deprecated(2019-04) - use subColumn()
        //  \deprecated(2019-04) - use subColumn()
        FOAM_DEPRECATED_FOR(2019-04, "subColumn() method")
        ConstMatrixBlock<mType> col
        (
            const label m,
            const label mStart,
            const label nStart
        ) const
        {
            return this->subColumn(nStart, mStart, m);
        }

        //- Deprecated(2019-04) - use subColumn()
        //  \deprecated(2019-04) - use subColumn()
        FOAM_DEPRECATED_FOR(2019-04, "subColumn() method")
        MatrixBlock<mType> col
        (
            const label m,
            const label mStart,
            const label nStart
        )
        {
            return this->subColumn(nStart, mStart, m);
        }

        //- Deleted(2019-04) - use subColumn()
        //  \deprecated(2019-04) - use subColumn()
        void col(const label m, const label rowStart) const = delete;

        //- Deleted(2019-04) - use subColumn()
        //  \deprecated(2019-04) - use subColumn()
        void col(const label m, const label rowStart) = delete;
};


// * * * * * * * * * * * * * * * IOstream Operator  * * * * * * * * * * * * * //

//- Read Matrix from Istream, discarding contents of existing Matrix.
template<class Form, class Type>
Istream& operator>>(Istream& is, Matrix<Form, Type>& mat)
{
    mat.readMatrix(is);
    return is;
}


//- Write Matrix to Ostream, as per Matrix::writeMatrix() with
//- default length, which is given by Detail::ListPolicy::short_length
template<class Form, class Type>
Ostream& operator<<(Ostream& os, const Matrix<Form, Type>& mat)
{
    return mat.writeMatrix(os, Detail::ListPolicy::short_length<Type>::value);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "MatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "Matrix.C"
    #include "MatrixIO.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
