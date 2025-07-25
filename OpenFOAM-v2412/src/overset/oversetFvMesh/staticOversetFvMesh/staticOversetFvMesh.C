/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify i
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

#include "staticOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(staticOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, staticOversetFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::staticOversetFvMesh::staticOversetFvMesh(const IOobject& io)
:
    staticFvMesh(io),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{
    oversetFvMeshBase::update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::staticOversetFvMesh::~staticOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::staticOversetFvMesh::update()
{
    if (!staticFvMesh::update())
    {
        return false;
    }

    oversetFvMeshBase::update();

    return true;
}


bool Foam::staticOversetFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok = staticFvMesh::writeObject(streamOpt, writeOnProc);
    ok = oversetFvMeshBase::writeObject(streamOpt, writeOnProc) && ok;
    return ok;
}


// ************************************************************************* //
