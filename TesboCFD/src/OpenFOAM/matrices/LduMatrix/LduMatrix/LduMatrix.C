/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "lduMatrix.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::LduMatrix(const lduMesh& mesh)
:
    lduMesh_(mesh)
{}


template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::LduMatrix(const LduMatrix& A)
:
    lduMesh_(A.lduMesh_)
{
    if (A.diagPtr_)
    {
        diagPtr_ = std::make_unique<Field<DType>>(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = std::make_unique<Field<LUType>>(*(A.upperPtr_));
    }

    if (A.lowerPtr_)
    {
        lowerPtr_ = std::make_unique<Field<LUType>>(*(A.lowerPtr_));
    }

    if (A.sourcePtr_)
    {
        sourcePtr_ = std::make_unique<Field<Type>>(*(A.sourcePtr_));
    }
}


template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::LduMatrix(LduMatrix&& A)
:
    lduMesh_(A.lduMesh_),
    diagPtr_(std::move(A.diagPtr_)),
    lowerPtr_(std::move(A.lowerPtr_)),
    upperPtr_(std::move(A.upperPtr_)),
    sourcePtr_(std::move(A.sourcePtr_))
{
    // Clear the old interfaces?
}


template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::LduMatrix(LduMatrix& A, bool reuse)
:
    lduMesh_(A.lduMesh_)
{
    if (reuse)
    {
        // Move assignment
        diagPtr_ = std::move(A.diagPtr_);
        upperPtr_ = std::move(A.upperPtr_);
        lowerPtr_ = std::move(A.lowerPtr_);
        sourcePtr_ = std::move(A.sourcePtr_);

        // Clear the old interfaces?
    }
    else
    {
        // Copy assignment
        if (A.diagPtr_)
        {
            diagPtr_ = std::make_unique<Field<DType>>(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = std::make_unique<Field<LUType>>(*(A.upperPtr_));
        }

        if (A.lowerPtr_)
        {
            lowerPtr_ = std::make_unique<Field<LUType>>(*(A.lowerPtr_));
        }

        if (A.sourcePtr_)
        {
            sourcePtr_ = std::make_unique<Field<Type>>(*(A.sourcePtr_));
        }
    }
}


template<class Type, class DType, class LUType>
Foam::LduMatrix<Type, DType, LUType>::LduMatrix
(
    const lduMesh& mesh,
    Istream& is
)
:
    lduMesh_(mesh),
    diagPtr_(new Field<DType>(is)),
    upperPtr_(new Field<LUType>(is)),
    lowerPtr_(new Field<LUType>(is)),
    sourcePtr_(new Field<Type>(is))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type, class DType, class LUType>
Foam::word Foam::LduMatrix<Type, DType, LUType>::matrixTypeName() const
{
    if (diagPtr_)
    {
        return
        (
            (!upperPtr_)
          ? (!lowerPtr_ ? "diagonal" : "diagonal-lower")
          : (!lowerPtr_ ? "symmetric" : "asymmetric")
        );
    }

    // is empty (or just wrong)
    return (!upperPtr_ && !lowerPtr_ ? "empty" : "ill-defined");
}


template<class Type, class DType, class LUType>
const Foam::Field<DType>& Foam::LduMatrix<Type, DType, LUType>::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorInFunction
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


template<class Type, class DType, class LUType>
Foam::Field<DType>& Foam::LduMatrix<Type, DType, LUType>::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ =
            std::make_unique<Field<DType>>(lduAddr().size(), Foam::zero{});
    }

    return *diagPtr_;
}


template<class Type, class DType, class LUType>
const Foam::Field<LUType>& Foam::LduMatrix<Type, DType, LUType>::upper() const
{
    if (upperPtr_)
    {
        return *upperPtr_;
    }
    else
    {
        if (!lowerPtr_)
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }

        return *lowerPtr_;
    }
}


template<class Type, class DType, class LUType>
Foam::Field<LUType>& Foam::LduMatrix<Type, DType, LUType>::upper()
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = std::make_unique<Field<LUType>>(*lowerPtr_);
        }
        else
        {
            upperPtr_ =
                std::make_unique<Field<LUType>>
                (
                    lduAddr().lowerAddr().size(),
                    Foam::zero{}
                );
        }
    }

    return *upperPtr_;
}


template<class Type, class DType, class LUType>
const Foam::Field<LUType>& Foam::LduMatrix<Type, DType, LUType>::lower() const
{
    if (lowerPtr_)
    {
        return *lowerPtr_;
    }
    else
    {
        if (!upperPtr_)
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }

        return *upperPtr_;
    }
}


template<class Type, class DType, class LUType>
Foam::Field<LUType>& Foam::LduMatrix<Type, DType, LUType>::lower()
{
    if (!lowerPtr_)
    {
        if (upperPtr_)
        {
            lowerPtr_ = std::make_unique<Field<LUType>>(*upperPtr_);
        }
        else
        {
            lowerPtr_ = std::make_unique<Field<LUType>>
            (
                lduAddr().lowerAddr().size(),
                Foam::zero{}
            );
        }
    }

    return *lowerPtr_;
}


template<class Type, class DType, class LUType>
const Foam::Field<Type>& Foam::LduMatrix<Type, DType, LUType>::source() const
{
    if (!sourcePtr_)
    {
        FatalErrorInFunction
            << "sourcePtr_ unallocated"
            << abort(FatalError);
    }

    return *sourcePtr_;
}


template<class Type, class DType, class LUType>
Foam::Field<Type>& Foam::LduMatrix<Type, DType, LUType>::source()
{
    if (!sourcePtr_)
    {
        sourcePtr_ =
            std::make_unique<Field<Type>>(lduAddr().size(), Foam::zero{});
    }

    return *sourcePtr_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// template<class Type, class DType, class LUType>
// Foam::Ostream& Foam::operator<<
// (
//     Ostream& os,
//     const InfoProxy<Type, DType, LUType>& iproxy
// )
// {
//     const auto& mat = *iproxy;
//
//     ...
//
//     os.check(FUNCTION_NAME);
//     return os;
// }


template<class Type, class DType, class LUType>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const LduMatrix<Type, DType, LUType>& mat
)
{
    if (mat.hasDiag())
    {
        os  << "Diagonal = " << mat.diag() << nl << nl;
    }

    if (mat.hasUpper())
    {
        os  << "Upper triangle = " << mat.upper() << nl << nl;
    }

    if (mat.hasLower())
    {
        os  << "Lower triangle = " << mat.lower() << nl << nl;
    }

    if (mat.hasSource())
    {
        os  << "Source = " << mat.source() << nl << nl;
    }

    os.check(FUNCTION_NAME);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "LduMatrixOperations.C"
#include "LduMatrixATmul.C"
#include "LduMatrixUpdateMatrixInterfaces.C"
#include "LduMatrixPreconditioner.C"
#include "LduMatrixSmoother.C"
#include "LduMatrixSolver.C"

// ************************************************************************* //
