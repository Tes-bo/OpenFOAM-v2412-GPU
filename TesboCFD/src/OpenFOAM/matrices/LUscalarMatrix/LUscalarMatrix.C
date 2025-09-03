/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2024 OpenCFD Ltd.
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

#include "LUscalarMatrix.H"
#include "lduMatrix.H"
#include "procLduMatrix.H"
#include "procLduInterface.H"
#include "cyclicLduInterface.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(LUscalarMatrix, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::LUscalarMatrix::LUscalarMatrix() noexcept
:
    comm_(UPstream::worldComm)
{}


Foam::LUscalarMatrix::LUscalarMatrix(const scalarSquareMatrix& mat)
:
    scalarSquareMatrix(mat),
    comm_(UPstream::worldComm)
{
    LUDecompose(*this, pivotIndices_);
}


Foam::LUscalarMatrix::LUscalarMatrix
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
:
    comm_(ldum.mesh().comm())
{
    if (UPstream::parRun())
    {
        PtrList<procLduMatrix> lduMatrices
        (
            UPstream::master(comm_) ? UPstream::nProcs(comm_) : 1
        );

        lduMatrices.set
        (
            0,  // rank-local matrix (and/or master)
            new procLduMatrix
            (
                ldum,
                interfaceCoeffs,
                interfaces
            )
        );

        if (UPstream::master(comm_))
        {
            for (const int proci : UPstream::subProcs(comm_))
            {
                auto& mat = lduMatrices.emplace_set(proci);

                IPstream::recv(mat, proci, UPstream::msgType(), comm_);
            }

            convert(lduMatrices);
        }
        else
        {
            OPstream::send
            (
                lduMatrices[0],  // rank-local matrix
                UPstream::masterNo(),
                UPstream::msgType(),
                comm_
            );
        }
    }
    else
    {
        convert(ldum, interfaceCoeffs, interfaces);
    }


    if (debug && UPstream::master(comm_))
    {
        const label numRows = nRows();
        const label numCols = nCols();

        Pout<< "LUscalarMatrix : size:" << numRows << endl;
        for (label rowi = 0; rowi < numRows; ++rowi)
        {
            const scalar* row = operator[](rowi);

            Pout<< "cell:" << rowi << " diagCoeff:" << row[rowi] << nl;

            Pout<< "    connects to upper cells :";
            for (label coli = rowi+1; coli < numCols; ++coli)
            {
                if (mag(row[coli]) > SMALL)
                {
                    Pout<< ' ' << coli << " (coeff:" << row[coli] << ')';
                }
            }
            Pout<< nl;
            Pout<< "    connects to lower cells :";
            for (label coli = 0; coli < rowi; ++coli)
            {
                if (mag(row[coli]) > SMALL)
                {
                    Pout<< ' ' << coli << " (coeff:" << row[coli] << ')';
                }
            }
            Pout<< nl;
        }
        Pout<< endl;
    }

    if (UPstream::master(comm_))
    {
        LUDecompose(*this, pivotIndices_);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::LUscalarMatrix::convert
(
    const lduMatrix& ldum,
    const FieldField<Field, scalar>& interfaceCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
)
{
    // Resize and fill with zero
    scalarSquareMatrix::resize_nocopy(ldum.lduAddr().size());
    scalarSquareMatrix::operator=(Foam::zero{});

    const label* __restrict__ uPtr = ldum.lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = ldum.lduAddr().lowerAddr().begin();

    const scalar* __restrict__ diagPtr = ldum.diag().begin();
    const scalar* __restrict__ upperPtr = ldum.upper().begin();
    const scalar* __restrict__ lowerPtr = ldum.lower().begin();

    const label nCells = ldum.diag().size();
    const label nFaces = ldum.upper().size();

    for (label cell=0; cell<nCells; cell++)
    {
        operator[](cell)[cell] = diagPtr[cell];
    }

    for (label face=0; face<nFaces; face++)
    {
        label uCell = uPtr[face];
        label lCell = lPtr[face];

        operator[](uCell)[lCell] = lowerPtr[face];
        operator[](lCell)[uCell] = upperPtr[face];
    }

    forAll(interfaces, inti)
    {
        if (interfaces.set(inti))
        {
            const lduInterface& interface = interfaces[inti].interface();

            // Assume any interfaces are cyclic ones

            const label* __restrict__ lPtr = interface.faceCells().begin();

            const cyclicLduInterface& cycInterface =
                refCast<const cyclicLduInterface>(interface);
            label nbrInt = cycInterface.neighbPatchID();
            const label* __restrict__ uPtr =
                interfaces[nbrInt].interface().faceCells().begin();

            const scalar* __restrict__ nbrUpperLowerPtr =
                interfaceCoeffs[nbrInt].begin();

            label inFaces = interface.faceCells().size();

            for (label face=0; face<inFaces; face++)
            {
                label uCell = lPtr[face];
                label lCell = uPtr[face];

                operator[](uCell)[lCell] -= nbrUpperLowerPtr[face];
            }
        }
    }
}


void Foam::LUscalarMatrix::convert
(
    const PtrList<procLduMatrix>& lduMatrices
)
{
    procOffsets_.resize_nocopy(lduMatrices.size() + 1);

    {
        auto iter = procOffsets_.begin();

        label nCellsTotal = 0;
        *iter++ = nCellsTotal;

        for (const auto& mat : lduMatrices)
        {
            nCellsTotal += mat.size();
            *iter++ = nCellsTotal;
        }

        // Resize and fill with zero
        scalarSquareMatrix::resize_nocopy(nCellsTotal);
        scalarSquareMatrix::operator=(Foam::zero{});
    }


    forAll(lduMatrices, ldumi)
    {
        const procLduMatrix& lduMatrixi = lduMatrices[ldumi];
        label offset = procOffsets_[ldumi];

        const label* __restrict__ uPtr = lduMatrixi.upperAddr_.begin();
        const label* __restrict__ lPtr = lduMatrixi.lowerAddr_.begin();

        const scalar* __restrict__ diagPtr = lduMatrixi.diag_.begin();
        const scalar* __restrict__ upperPtr = lduMatrixi.upper_.begin();
        const scalar* __restrict__ lowerPtr = lduMatrixi.lower_.begin();

        const label nCells = lduMatrixi.size();
        const label nFaces = lduMatrixi.upper_.size();

        for (label cell=0; cell<nCells; cell++)
        {
            label globalCell = cell + offset;
            operator[](globalCell)[globalCell] = diagPtr[cell];
        }

        for (label face=0; face<nFaces; face++)
        {
            label uCell = uPtr[face] + offset;
            label lCell = lPtr[face] + offset;

            operator[](uCell)[lCell] = lowerPtr[face];
            operator[](lCell)[uCell] = upperPtr[face];
        }

        const PtrList<procLduInterface>& interfaces =
            lduMatrixi.interfaces_;

        forAll(interfaces, inti)
        {
            const procLduInterface& interface = interfaces[inti];

            if (interface.myProcNo_ == interface.neighbProcNo_)
            {
                const label* __restrict__ ulPtr = interface.faceCells_.begin();

                const scalar* __restrict__ upperLowerPtr =
                    interface.coeffs_.begin();

                label inFaces = interface.faceCells_.size()/2;

                for (label face=0; face<inFaces; face++)
                {
                    label uCell = ulPtr[face] + offset;
                    label lCell = ulPtr[face + inFaces] + offset;

                    operator[](uCell)[lCell] -= upperLowerPtr[face + inFaces];
                    operator[](lCell)[uCell] -= upperLowerPtr[face];
                }
            }
            else if (interface.myProcNo_ < interface.neighbProcNo_)
            {
                // Interface to neighbour proc. Find on neighbour proc the
                // corresponding interface. The problem is that there can
                // be multiple interfaces between two processors (from
                // processorCyclics) so also compare the communication tag

                const PtrList<procLduInterface>& neiInterfaces =
                    lduMatrices[interface.neighbProcNo_].interfaces_;

                label neiInterfacei = -1;

                forAll(neiInterfaces, ninti)
                {
                    if
                    (
                        (
                            neiInterfaces[ninti].neighbProcNo_
                         == interface.myProcNo_
                        )
                     && (neiInterfaces[ninti].tag_ ==  interface.tag_)
                    )
                    {
                        neiInterfacei = ninti;
                        break;
                    }
                }

                if (neiInterfacei == -1)
                {
                    FatalErrorInFunction << exit(FatalError);
                }

                const procLduInterface& neiInterface =
                    neiInterfaces[neiInterfacei];

                const label* __restrict__ uPtr = interface.faceCells_.begin();
                const label* __restrict__ lPtr =
                    neiInterface.faceCells_.begin();

                const scalar* __restrict__ upperPtr = interface.coeffs_.begin();
                const scalar* __restrict__ lowerPtr =
                    neiInterface.coeffs_.begin();

                label inFaces = interface.faceCells_.size();
                label neiOffset = procOffsets_[interface.neighbProcNo_];

                for (label face=0; face<inFaces; face++)
                {
                    label uCell = uPtr[face] + offset;
                    label lCell = lPtr[face] + neiOffset;

                    operator[](uCell)[lCell] -= lowerPtr[face];
                    operator[](lCell)[uCell] -= upperPtr[face];
                }
            }
        }
    }
}


void Foam::LUscalarMatrix::printDiagonalDominance() const
{
    for (label i=0; i<m(); i++)
    {
        scalar sum = 0.0;
        for (label j=0; j<m(); j++)
        {
            if (i != j)
            {
                sum += operator[](i)[j];
            }
        }
        Info<< mag(sum)/mag(operator[](i)[i]) << endl;
    }
}


void Foam::LUscalarMatrix::decompose(const scalarSquareMatrix& mat)
{
    scalarSquareMatrix::operator=(mat);
    LUDecompose(*this, pivotIndices_);
}


void Foam::LUscalarMatrix::inv(scalarSquareMatrix& M) const
{
    scalarField source(m());

    for (label j=0; j<m(); j++)
    {
        source = Zero;
        source[j] = 1;
        LUBacksubstitute(*this, pivotIndices_, source);
        for (label i=0; i<m(); i++)
        {
            M(i, j) = source[i];
        }
    }
}


// ************************************************************************* //
