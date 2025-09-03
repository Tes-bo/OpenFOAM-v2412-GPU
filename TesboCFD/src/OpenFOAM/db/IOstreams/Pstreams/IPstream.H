/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2021-2024 OpenCFD Ltd.
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
    Foam::IPstream

Description
    Input inter-processor communications stream.

SourceFiles
    IPstreams.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef Foam_IPstream_H
#define Foam_IPstream_H

#include "UIPstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IPstream Declaration
\*---------------------------------------------------------------------------*/

class IPstream
:
    public Pstream,
    public UIPstream
{
public:

    // Constructors

        //- Construct given process index to read from
        //- and optional buffer size, read format
        IPstream
        (
            const UPstream::commsTypes commsType,
            const int fromProcNo,
            const label bufSize = 0,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    // Static Functions

        //- Receive and deserialize a value.
        //- Uses \c operator>> for de-serialization
        template<class Type>
        static void recv
        (
            Type& value,
            const int fromProcNo,
            const int tag = UPstream::msgType(),
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        {
            IPstream is
            (
                UPstream::commsTypes::scheduled,  // ie, MPI_Recv()
                fromProcNo,
                0,  // bufSize
                tag,
                comm,
                fmt
            );
            is >> value;
        }
};


/*---------------------------------------------------------------------------*\
                          Class IPBstream Declaration
\*---------------------------------------------------------------------------*/

//- Input inter-processor communications stream
//- using MPI broadcast.
class IPBstream
:
    public Pstream,
    public UIPBstream
{
public:

    // Constructors

        //- Construct for broadcast root, optional buffer size, read format
        IPBstream
        (
            const UPstream::commsTypes,     //!< ignored
            const int rootProcNo,           //!< normally UPstream::masterNo()
            const label bufSize = 0,
            const int tag = UPstream::msgType(),  //!< ignored
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct for broadcast root and communicator,
        //- with optional read format
        IPBstream
        (
            const int rootProcNo,           //!< normally UPstream::masterNo()
            const label comm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );

        //- Construct with optional communicator and read format.
        //- Uses UPstream::masterNo() root
        explicit IPBstream
        (
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        );


    // Static Functions

        //- Receive (from broadcast, root == UPstream::masterNo())
        //- and deserialize a value.
        //- Uses \c operator>> for de-serialization
        template<class Type>
        static void recv
        (
            Type& value,
            const label comm = UPstream::worldComm,
            IOstreamOption::streamFormat fmt = IOstreamOption::BINARY
        )
        {
            IPBstream is(comm, fmt);
            is >> value;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
