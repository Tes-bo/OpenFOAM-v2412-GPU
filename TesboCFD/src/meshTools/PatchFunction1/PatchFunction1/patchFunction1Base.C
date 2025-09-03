/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2021 OpenCFD Ltd.
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

#include "patchFunction1Base.H"
#include "polyBoundaryMesh.H"
#include "polyMesh.H"
#include "polyPatch.H"
#include "objectRegistry.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::patchFunction1Base::patchFunction1Base
(
    const polyPatch& pp,
    const word& entryName,
    const bool faceValues
)
:
    refCount(),
    name_(entryName),
    patch_(pp),
    faceValues_(faceValues)
{}


Foam::patchFunction1Base::patchFunction1Base
(
    const polyPatch& pp,
    const word& entryName,
    const dictionary& dict,
    const bool faceValues
)
:
    refCount(),
    name_(entryName),
    patch_(pp),
    faceValues_(faceValues)
{}


Foam::patchFunction1Base::patchFunction1Base(const patchFunction1Base& rhs)
:
    patchFunction1Base(rhs, rhs.patch())
{}


Foam::patchFunction1Base::patchFunction1Base
(
    const patchFunction1Base& rhs,
    const polyPatch& pp
)
:
    refCount(),
    name_(rhs.name_),
    patch_(pp),
    faceValues_(rhs.faceValues_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry* Foam::patchFunction1Base::whichDb() const
{
    return &(patch_.boundaryMesh().mesh());  // mesh registry
}


const Foam::objectRegistry& Foam::patchFunction1Base::obr() const
{
    return patch_.boundaryMesh().mesh();  // mesh registry
}


const Foam::Time& Foam::patchFunction1Base::time() const
{
    return patch_.boundaryMesh().mesh().time();
}


void Foam::patchFunction1Base::userTimeToTime(const Time& t)
{}


// ************************************************************************* //
