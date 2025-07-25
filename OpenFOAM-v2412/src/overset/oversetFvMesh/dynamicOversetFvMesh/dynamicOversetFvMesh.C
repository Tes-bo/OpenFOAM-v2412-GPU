/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022,2024 OpenCFD Ltd.
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

#include "dynamicOversetFvMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dynamicOversetFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, dynamicOversetFvMesh, IOobject);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::dynamicOversetFvMesh
(
    const IOobject& io,
    const bool doInit
)
:
    dynamicMotionSolverListFvMesh(io, doInit),
    oversetFvMeshBase(static_cast<const fvMesh&>(*this))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dynamicOversetFvMesh::~dynamicOversetFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dynamicOversetFvMesh::update()
{
    if (!dynamicMotionSolverListFvMesh::update())
    {
        return false;
    }

    //Note: too late to do oversetFvMeshBase::clearOut() to get it
    //      consistent with any new cell-cell stencil since 
    //      dynamicMotionSolverListFvMesh already triggers
    //      meshObject::movePoints on cellCellStencilObject
    oversetFvMeshBase::update();

    return true;
}


bool Foam::dynamicOversetFvMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    bool ok =
        dynamicMotionSolverListFvMesh::writeObject(streamOpt, writeOnProc);
    ok = oversetFvMeshBase::writeObject(streamOpt, writeOnProc) && ok;
    return ok;
}


// ************************************************************************* //
