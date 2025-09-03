/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2022 OpenCFD Ltd.
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

#include "MeshedSurfaceIOAllocator.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints,
    const IOobject& ioFaces
)
:
    points_(ioPoints),
    faces_(ioFaces)
{}


Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints, const pointField& points,
    const IOobject& ioFaces,  const faceList& faces
)
:
    points_(ioPoints, points),
    faces_(ioFaces, faces)
{}


Foam::Detail::MeshedSurfaceIOAllocator::MeshedSurfaceIOAllocator
(
    const IOobject& ioPoints, pointField&& points,
    const IOobject& ioFaces,  faceList&& faces
)
:
    points_(ioPoints, std::move(points)),
    faces_(ioFaces, std::move(faces))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::Detail::MeshedSurfaceIOAllocator::setInstance
(
    const fileName& inst
)
{
    points_.instance() = inst;
    faces_.instance()  = inst;
}


void Foam::Detail::MeshedSurfaceIOAllocator::setWriteOption
(
    IOobjectOption::writeOption wOpt
)
{
    points_.writeOpt(wOpt);
    faces_.writeOpt(wOpt);
}


void Foam::Detail::MeshedSurfaceIOAllocator::clear()
{
    points_.clear();
    faces_.clear();
}


bool Foam::Detail::MeshedSurfaceIOAllocator::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    return
    (
        points_.writeObject(streamOpt, writeOnProc)
     && faces_.writeObject(streamOpt, writeOnProc)
    );
}


// ************************************************************************* //
