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

#include "LUscalarMatrix.H"
#include "SubList.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::LUscalarMatrix::solve
(
    List<Type>& x,
    const UList<Type>& source
) const
{
    // If x and source are different initialize x = source
    if (&x != &source)
    {
        x = source;
    }

    const auto tag = UPstream::msgType();

    if (UPstream::parRun())
    {
        List<Type> allx;  // scratch space (on master)

        const label startOfRequests = UPstream::nRequests();

        // Like globalIndex::gather()
        if (UPstream::master(comm_))
        {
            allx.resize(m());

            SubList<Type>(allx, x.size()) = x;

            for (const int proci : UPstream::subProcs(comm_))
            {
                SubList<Type> procSlot
                (
                    allx,
                    procOffsets_[proci+1]-procOffsets_[proci],
                    procOffsets_[proci]
                );

                if (procSlot.empty())
                {
                    // Nothing to do
                }
                else if (is_contiguous<Type>::value)
                {
                    UIPstream::read
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlot.data_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm_
                    );
                }
                else
                {
                    IPstream::recv(procSlot, proci, tag, comm_);
                }
            }
        }
        else
        {
            if (x.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UOPstream::write
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    x.cdata_bytes(),
                    x.size_bytes(),
                    tag,
                    comm_
                );
            }
            else
            {
                OPstream::send(x, UPstream::masterNo(), tag, comm_);
            }
        }

        UPstream::waitRequests(startOfRequests);

        // LUBacksubstitute and then like globalIndex::scatter()
        if (UPstream::master(comm_))
        {
            LUBacksubstitute(*this, pivotIndices_, allx);

            x = SubList<Type>(allx, x.size());

            for (const int proci : UPstream::subProcs(comm_))
            {
                SubList<Type> procSlot
                (
                    allx,
                    procOffsets_[proci+1]-procOffsets_[proci],
                    procOffsets_[proci]
                );

                if (procSlot.empty())
                {
                    // Nothing to do
                }
                else if (is_contiguous<Type>::value)
                {
                    UOPstream::write
                    (
                        UPstream::commsTypes::nonBlocking,
                        proci,
                        procSlot.cdata_bytes(),
                        procSlot.size_bytes(),
                        tag,
                        comm_
                    );
                }
                else
                {
                    OPstream::send(procSlot, proci, tag, comm_);
                }
            }
        }
        else
        {
            if (x.empty())
            {
                // Nothing to do
            }
            else if (is_contiguous<Type>::value)
            {
                UIPstream::read
                (
                    UPstream::commsTypes::nonBlocking,
                    UPstream::masterNo(),
                    x.data_bytes(),
                    x.size_bytes(),
                    tag,
                    comm_
                );
            }
            else
            {
                IPstream::recv(x, UPstream::masterNo(), tag, comm_);
            }
        }

        UPstream::waitRequests(startOfRequests);
    }
    else
    {
        LUBacksubstitute(*this, pivotIndices_, x);
    }
}


template<class Type>
Foam::tmp<Foam::Field<Type>> Foam::LUscalarMatrix::solve
(
    const UList<Type>& source
) const
{
    auto tx = tmp<Field<Type>>::New(m());

    solve(tx.ref(), source);

    return tx;
}


// ************************************************************************* //
