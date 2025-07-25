/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2024 OpenCFD Ltd.
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
#include "Switch.H"
#include "objectRegistry.H"
#include "scalarIOField.H"
#include "Time.H"
#include "meshState.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(lduMatrix, 1);
}


const Foam::scalar Foam::lduMatrix::defaultTolerance = 1e-6;

const Foam::Enum
<
    Foam::lduMatrix::normTypes
>
Foam::lduMatrix::normTypesNames_
({
    { normTypes::NO_NORM, "none" },
    { normTypes::DEFAULT_NORM, "default" },
    { normTypes::L1_SCALED_NORM, "L1_scaled" },
});


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::lduMatrix::lduMatrix(const lduMesh& mesh)
:
    lduMesh_(mesh)
{}


Foam::lduMatrix::lduMatrix(const lduMatrix& A)
:
    lduMesh_(A.lduMesh_)
{
    if (A.diagPtr_)
    {
        diagPtr_ = std::make_unique<scalarField>(*(A.diagPtr_));
    }

    if (A.upperPtr_)
    {
        upperPtr_ = std::make_unique<scalarField>(*(A.upperPtr_));
    }

    if (A.lowerPtr_)
    {
        lowerPtr_ = std::make_unique<scalarField>(*(A.lowerPtr_));
    }

    if (A.lowerCSRPtr_)
    {
        lowerCSRPtr_ = std::make_unique<scalarField>(*(A.lowerCSRPtr_));
    }

    // No need to copy work
}


Foam::lduMatrix::lduMatrix(lduMatrix&& A)
:
    lduMesh_(A.lduMesh_),
    diagPtr_(std::move(A.diagPtr_)),
    lowerPtr_(std::move(A.lowerPtr_)),
    upperPtr_(std::move(A.upperPtr_)),
    lowerCSRPtr_(std::move(A.lowerCSRPtr_)),
    workPtr_(std::move(A.workPtr_))
{}


Foam::lduMatrix::lduMatrix(lduMatrix& A, bool reuse)
:
    lduMesh_(A.lduMesh_)
{
    if (reuse)
    {
        diagPtr_ = std::move(A.diagPtr_);
        upperPtr_ = std::move(A.upperPtr_);
        lowerPtr_ = std::move(A.lowerPtr_);
        lowerCSRPtr_ = std::move(A.lowerCSRPtr_);
        workPtr_ = std::move(A.workPtr_);
    }
    else
    {
        // Copy assignment
        if (A.diagPtr_)
        {
            diagPtr_ = std::make_unique<scalarField>(*(A.diagPtr_));
        }

        if (A.upperPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*(A.upperPtr_));
        }

        if (A.lowerPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*(A.lowerPtr_));
        }

        // Note: no real need to keep lowerCSR except we use it (hasLowerCSR())
        //       to trigger certain actions
        if (A.lowerCSRPtr_)
        {
            lowerCSRPtr_ = std::make_unique<scalarField>(*(A.lowerCSRPtr_));
        }
    }
}


Foam::lduMatrix::lduMatrix(const lduMesh& mesh, Istream& is)
:
    lduMesh_(mesh)
{
    bool withLower, withDiag, withUpper, withLowerCSR;

    is >> withLower >> withDiag >> withUpper >> withLowerCSR;

    if (withLower)
    {
        lowerPtr_ = std::make_unique<scalarField>(is);
    }
    if (withDiag)
    {
        diagPtr_ = std::make_unique<scalarField>(is);
    }
    if (withUpper)
    {
        upperPtr_ = std::make_unique<scalarField>(is);
    }
    if (withLowerCSR)
    {
        // Check if it has been sent over or we need to construct it. We need
        // to set the lowerCSRPtr since it determines behaviour.
        if (withLower)
        {
            // Construct from lower
            (void)lowerCSR();
        }
        else
        {
            lowerCSRPtr_ = std::make_unique<scalarField>(is);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::word Foam::lduMatrix::matrixTypeName() const
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


const Foam::scalarField& Foam::lduMatrix::diag() const
{
    if (!diagPtr_)
    {
        FatalErrorInFunction
            << "diagPtr_ unallocated"
            << abort(FatalError);
    }

    return *diagPtr_;
}


Foam::scalarField& Foam::lduMatrix::diag()
{
    if (!diagPtr_)
    {
        diagPtr_ =
            std::make_unique<scalarField>(lduAddr().size(), Foam::zero{});
    }

    return *diagPtr_;
}


Foam::scalarField& Foam::lduMatrix::diag(label size)
{
    if (!diagPtr_)
    {
        // if (size < 0)
        // {
        //     size = lduAddr().size();
        // }
        diagPtr_ = std::make_unique<scalarField>(size, Foam::zero{});
    }

    return *diagPtr_;
}


const Foam::scalarField& Foam::lduMatrix::upper() const
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


Foam::scalarField& Foam::lduMatrix::upper()
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*lowerPtr_);
        }
        else
        {
            // no lowerPtr so any lowerCSR was constructed from upper
            lowerCSRPtr_.reset(nullptr);

            upperPtr_ =
                std::make_unique<scalarField>
                (
                    lduAddr().lowerAddr().size(),
                    Foam::zero{}
                );
        }
    }

    return *upperPtr_;
}


Foam::scalarField& Foam::lduMatrix::upper(label nCoeffs)
{
    if (!upperPtr_)
    {
        if (lowerPtr_)
        {
            upperPtr_ = std::make_unique<scalarField>(*lowerPtr_);
        }
        else
        {
            // no lowerPtr so any lowerCSR was constructed from upper
            lowerCSRPtr_.reset(nullptr);

            // if (nCoeffs < 0)
            // {
            //     nCoeffs = lduAddr().lowerAddr().size();
            // }
            upperPtr_ = std::make_unique<scalarField>(nCoeffs, Foam::zero{});
        }
    }

    return *upperPtr_;
}


const Foam::scalarField& Foam::lduMatrix::lower() const
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


Foam::scalarField& Foam::lduMatrix::lower()
{
    if (!lowerPtr_)
    {
        lowerCSRPtr_.reset(nullptr);

        if (upperPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*upperPtr_);
        }
        else
        {
            lowerPtr_ =
                std::make_unique<scalarField>
                (
                    lduAddr().lowerAddr().size(),
                    Foam::zero{}
                );
        }
    }

    return *lowerPtr_;
}


Foam::scalarField& Foam::lduMatrix::lower(label nCoeffs)
{
    if (!lowerPtr_)
    {
        lowerCSRPtr_.reset(nullptr);

        if (upperPtr_)
        {
            lowerPtr_ = std::make_unique<scalarField>(*upperPtr_);
        }
        else
        {
            // if (nCoeffs < 0)
            // {
            //     nCoeffs = lduAddr().lowerAddr().size();
            // }
            lowerPtr_ =
                std::make_unique<scalarField>(nCoeffs, Foam::zero{});
        }
    }

    return *lowerPtr_;
}


const Foam::scalarField& Foam::lduMatrix::lowerCSR() const
{
    if (!lowerCSRPtr_)
    {
        const label nLower = lduAddr().losortAddr().size();

        lowerCSRPtr_ = std::make_unique<scalarField>(nLower);

        if (lowerPtr_)
        {
            lduAddr().map(*lowerPtr_, *lowerCSRPtr_);
        }
        else if (upperPtr_)
        {
            lduAddr().map(*upperPtr_, *lowerCSRPtr_);
        }
        else
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }
    }

    return *lowerCSRPtr_;
}


Foam::scalarField& Foam::lduMatrix::lowerCSR()
{
    if (!lowerCSRPtr_)
    {
        const label nLower = lduAddr().losortAddr().size();

        lowerCSRPtr_ = std::make_unique<scalarField>(nLower);

        if (lowerPtr_)
        {
            lduAddr().map(*lowerPtr_, *lowerCSRPtr_);
        }
        else if (upperPtr_)
        {
            lduAddr().map(*upperPtr_, *lowerCSRPtr_);
        }
        else
        {
            FatalErrorInFunction
                << "lowerPtr_ and upperPtr_ unallocated"
                << abort(FatalError);
        }
    }

    return *lowerCSRPtr_;
}


Foam::solveScalarField& Foam::lduMatrix::work(label size) const
{
    if (!workPtr_ || workPtr_->size() != size)
    {
        workPtr_ = std::make_unique<solveScalarField>(size, Foam::zero{});
    }

    return *workPtr_;
}


const Foam::solveScalarField& Foam::lduMatrix::work() const
{
    if (!workPtr_)
    {
        FatalErrorInFunction
            << "workPtr_ unallocated"
            << abort(FatalError);
    }

    return *workPtr_;
}


void Foam::lduMatrix::setResidualField
(
    const scalarField& residual,
    const word& fieldName,
    const bool initial
) const
{
    if (!mesh().hasDb())
    {
        return;
    }

    scalarIOField* residualPtr =
        mesh().thisDb().getObjectPtr<scalarIOField>
        (
            initial
          ? IOobject::scopedName("initialResidual", fieldName)
          : IOobject::scopedName("residual", fieldName)
        );

    if (residualPtr)
    {
        const auto* dataPtr = mesh().thisDb().findObject<meshState>("data");

        if (dataPtr)
        {
            if (initial && dataPtr->isFirstIteration())
            {
                *residualPtr = residual;
                DebugInfo
                    << "Setting residual field for first solver iteration "
                    << "for solver field: " << fieldName << endl;
            }
        }
        else
        {
            *residualPtr = residual;
            DebugInfo
                << "Setting residual field for solver field "
                << fieldName << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const lduMatrix& mat)
{
    os  << mat.hasLower() << token::SPACE
        << mat.hasDiag() << token::SPACE
        << mat.hasUpper() << token::SPACE
        << mat.hasLowerCSR() << token::SPACE;

    if (mat.hasLower())
    {
        os  << mat.lower();
    }

    if (mat.hasDiag())
    {
        os  << mat.diag();
    }

    if (mat.hasUpper())
    {
        os  << mat.upper();
    }

    if (mat.hasLowerCSR() && !mat.hasLower())
    {
        // Only send over if can not be reconstructed locally
        os  << mat.lowerCSR();
    }

    os.check(FUNCTION_NAME);

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const InfoProxy<lduMatrix>& iproxy
)
{
    const auto& mat = *iproxy;

    os  << "Lower:" << Switch::name(mat.hasLower())
        << " Diag:" << Switch::name(mat.hasDiag())
        << " Upper:" << Switch::name(mat.hasUpper())
        << " lowerCSR:" << Switch::name(mat.hasLowerCSR())
        << endl;

    if (mat.hasLower())
    {
        os  << "lower:" << mat.lower().size() << endl;
    }
    if (mat.hasDiag())
    {
        os  << "diag :" << mat.diag().size() << endl;
    }
    if (mat.hasUpper())
    {
        os  << "upper:" << mat.upper().size() << endl;
    }


    //if (hasLower)
    //{
    //    os  << "lower contents:" << endl;
    //    forAll(mat.lower(), i)
    //    {
    //        os  << "i:" << i << "\t" << mat.lower()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasDiag)
    //{
    //    os  << "diag contents:" << endl;
    //    forAll(mat.diag(), i)
    //    {
    //        os  << "i:" << i << "\t" << mat.diag()[i] << endl;
    //    }
    //    os  << endl;
    //}
    //if (hasUpper)
    //{
    //    os  << "upper contents:" << endl;
    //    forAll(mat.upper(), i)
    //    {
    //        os  << "i:" << i << "\t" << mat.upper()[i] << endl;
    //    }
    //    os  << endl;
    //}

    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
