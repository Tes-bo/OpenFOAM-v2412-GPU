/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023-2024 OpenCFD Ltd.
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

#include "faMeshesRegistry.H"
#include "faMesh.H"
#include "polyMesh.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faMeshesRegistry, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::faMeshesRegistry::faMeshesRegistry(const polyMesh& mesh)
:
    MeshObject_type(mesh),
    objects_
    (
        IOobject
        (
            faMesh::prefix(),
            mesh.time().timeName(),
            mesh.thisDb(),
            IOobjectOption::NO_READ,
            IOobjectOption::AUTO_WRITE,
            IOobjectOption::REGISTER
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faMeshesRegistry::movePoints()
{
    for (faMesh& m : objects_.sorted<faMesh>())
    {
        m.movePoints();
    }

    return true;
}


void Foam::faMeshesRegistry::updateMesh(const mapPolyMesh& mpm)
{
    for (faMesh& m : objects_.sorted<faMesh>())
    {
        m.updateMesh(mpm);
    }
}


bool Foam::faMeshesRegistry::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    // for (const faMesh& m : objects_.csorted<faMesh>())
    // {
    //     m.write(writeOnProc);
    // }
    //
    // return true;

    return objects_.writeObject(streamOpt, writeOnProc);
}


// ************************************************************************* //
