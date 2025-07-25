/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022-2024 OpenCFD Ltd.
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

#include "OPstream.H"
#include "IPstream.H"
#include "contiguous.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
void Foam::Pstream::broadcast(Type& value, const label comm)
{
    if (is_contiguous<Type>::value)
    {
        // Note: contains parallel guard internally
        UPstream::broadcast
        (
            reinterpret_cast<char*>(&value),
            sizeof(Type),
            comm
        );
    }
    else if (UPstream::is_parallel(comm))
    {
        if (UPstream::master(comm))
        {
            OPBstream os(comm);
            os << value;
        }
        else  // UPstream::is_subrank(comm)
        {
            IPBstream is(comm);
            is >> value;
        }
    }
}


template<class Type, class... Args>
void Foam::Pstream::broadcasts(const label comm, Type& arg1, Args&&... args)
{
    if (UPstream::is_parallel(comm))
    {
        if (UPstream::master(comm))
        {
            OPBstream os(comm);
            Detail::outputLoop(os, arg1, std::forward<Args>(args)...);
        }
        else  // UPstream::is_subrank(comm)
        {
            IPBstream is(comm);
            Detail::inputLoop(is, arg1, std::forward<Args>(args)...);
        }
    }
}


template<class ListType>
void Foam::Pstream::broadcastList(ListType& list, const label comm)
{
    if (is_contiguous<typename ListType::value_type>::value)
    {
        // List data are contiguous
        // 1. broadcast the size
        // 2. resize for receiver list
        // 3. broadcast contiguous contents

        if (UPstream::is_parallel(comm))
        {
            label len(list.size());

            UPstream::broadcast
            (
                reinterpret_cast<char*>(&len),
                sizeof(label),
                comm
            );

            if (UPstream::is_subrank(comm))
            {
                list.resize_nocopy(len);
            }

            if (len)
            {
                UPstream::broadcast
                (
                    list.data_bytes(),
                    list.size_bytes(),
                    comm
                );
            }
        }
    }
    else if (UPstream::is_parallel(comm))
    {
        // List data are non-contiguous - serialize/de-serialize

        if (UPstream::master(comm))
        {
            OPBstream os(comm);
            os << list;
        }
        else  // UPstream::is_subrank(comm)
        {
            IPBstream is(comm);
            is >> list;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Convenience wrappers - defined after all specialisations are known

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a broadcasted value (uses a copy internally)
template<class Type>
Type returnBroadcast
(
    const Type& value,
    const label comm = UPstream::worldComm
)
{
    Type work(value);
    Pstream::broadcast(work, comm);
    return work;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// ************************************************************************* //
