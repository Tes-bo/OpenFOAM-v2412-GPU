/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015 OpenFOAM Foundation
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

#include "IOmapDistributePolyMesh.H"

/* * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * */

namespace Foam
{
    defineTypeNameAndDebug(IOmapDistributePolyMesh, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::IOmapDistributePolyMesh::readIOcontents()
{
    if (isReadRequired() || (isReadOptional() && headerOk()))
    {
        readStream(typeName) >> *this;
        close();
        return true;
    }

    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOmapDistributePolyMesh::IOmapDistributePolyMesh(const IOobject& io)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistributePolyMesh>();

    readIOcontents();
}


Foam::IOmapDistributePolyMesh::IOmapDistributePolyMesh
(
    const IOobject& io,
    const mapDistributePolyMesh& map
)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistributePolyMesh>();

    if (!readIOcontents())
    {
        mapDistributePolyMesh::operator=(map);
    }
}


Foam::IOmapDistributePolyMesh::IOmapDistributePolyMesh
(
    const IOobject& io,
    mapDistributePolyMesh&& map
)
:
    regIOobject(io)
{
    // Warn for MUST_READ_IF_MODIFIED
    warnNoRereading<IOmapDistributePolyMesh>();

    mapDistributePolyMesh::transfer(map);

    readIOcontents();
}


Foam::IOmapDistributePolyMeshRef::IOmapDistributePolyMeshRef
(
    const IOobject& io,
    const mapDistributePolyMesh& map
)
:
    regIOobject(io),
    contentRef_(map)  // cref
{}


// Not sure if we need this yet...
//
/// Foam::IOmapDistributePolyMeshRef::IOmapDistributePolyMeshRef
/// (
///     const IOobject& io,
///     mapDistributePolyMesh& map
/// )
/// :
///     regIOobject(io),
///     contentRef_()
/// {
///     contentRef_.ref(map);  // writable reference
/// }


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::IOmapDistributePolyMesh::readData(Istream& is)
{
    is >> *this;
    return is.good();
}


bool Foam::IOmapDistributePolyMesh::writeData(Ostream& os) const
{
    os << *this;
    return os.good();
}


bool Foam::IOmapDistributePolyMeshRef::readData(Istream& is)
{
    is >> contentRef_.ref();
    return is.good();
}


bool Foam::IOmapDistributePolyMeshRef::writeData(Ostream& os) const
{
    os << contentRef_.cref();
    return os.good();
}


// ************************************************************************* //
