/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2018 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "threadedCollatedOFstream.H"
#include "decomposedBlockData.H"
#include "OFstreamCollator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::threadedCollatedOFstream::threadedCollatedOFstream
(
    OFstreamCollator& writer,
    IOstreamOption::atomicType atomic,
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool useThread
)
:
    OStringStream(streamOpt),
    writer_(writer),
    pathName_(pathName),
    atomic_(atomic),
    compression_(streamOpt.compression()),
    useThread_(useThread),
    headerEntries_()
{}


Foam::threadedCollatedOFstream::threadedCollatedOFstream
(
    OFstreamCollator& writer,
    const fileName& pathName,
    IOstreamOption streamOpt,
    const bool useThread
)
:
    threadedCollatedOFstream
    (
        writer,
        IOstreamOption::NON_ATOMIC,
        pathName,
        streamOpt,
        useThread
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::threadedCollatedOFstream::~threadedCollatedOFstream()
{
    writer_.write
    (
        decomposedBlockData::typeName,
        pathName_,
        str(),
        IOstreamOption(IOstreamOption::BINARY, version(), compression_),
        atomic_,
        IOstreamOption::NO_APPEND,
        useThread_,
        headerEntries_
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::threadedCollatedOFstream::setHeaderEntries(const dictionary& dict)
{
    headerEntries_ = dict;
}


// ************************************************************************* //
