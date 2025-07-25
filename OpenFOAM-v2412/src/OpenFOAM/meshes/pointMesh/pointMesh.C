/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "pointMesh.H"
#include "globalMeshData.H"
#include "pointMeshMapper.H"
#include "pointFields.H"
#include "MapGeometricFields.H"
#include "MapPointField.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pointMesh, 0);
}

Foam::word Foam::pointMesh::meshSubDir = "pointMesh";


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointMesh::mapFields(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "void pointMesh::mapFields(const mapPolyMesh&): "
            << "Mapping all registered pointFields."
            << endl;
    }
    // Create a mapper
    const pointMeshMapper m(*this, mpm);

    MapGeometricFields<scalar, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields<vector, pointPatchField, pointMeshMapper, pointMesh>(m);
    MapGeometricFields
    <
        sphericalTensor,
        pointPatchField,
        pointMeshMapper,
        pointMesh
    >(m);
    MapGeometricFields<symmTensor, pointPatchField, pointMeshMapper, pointMesh>
    (m);
    MapGeometricFields<tensor, pointPatchField, pointMeshMapper, pointMesh>(m);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointMesh::pointMesh(const polyMesh& pMesh)
:
    MeshObject_type(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(*this, pMesh.boundaryMesh())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from polyMesh " << pMesh.name()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


Foam::pointMesh::pointMesh(const polyMesh& pMesh, const IOobject& io)
:
    MeshObject<polyMesh, Foam::UpdateableMeshObject, pointMesh>(pMesh),
    GeoMesh<polyMesh>(pMesh),
    boundary_(io, *this, pMesh.boundaryMesh())
{
    if (debug)
    {
        Pout<< "pointMesh::pointMesh(const polyMesh&): "
            << "Constructing from IO " << io.objectRelPath()
            << endl;
    }

    // Calculate the geometry for the patches (transformation tensors etc.)
    boundary_.calcGeometry();
}


Foam::pointMesh::pointMesh
(
    const polyMesh& pMesh,
    const IOobjectOption::readOption rOpt
)
:
    pointMesh
    (
        pMesh,
        IOobject
        (
            pMesh.name(),                // polyMesh region
            pMesh.facesInstance(),       // polyMesh topology instance
            pMesh.time(),
            rOpt,
            Foam::IOobject::NO_WRITE
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pointMesh::setInstance
(
    const fileName& inst,
    const IOobjectOption::writeOption wOpt
)
{
    if (debug)
    {
        Pout<< "pointMesh::setInstance(): "
            << "Setting instance to " << inst << endl;
    }
    this->writeOpt(wOpt);
    this->instance() = inst;

    boundary_.writeOpt(wOpt);
    boundary_.instance() = inst;
}


bool Foam::pointMesh::movePoints()
{
    if (debug)
    {
        Pout<< "pointMesh::movePoints(): "
            << "Moving points." << endl;
    }

    boundary_.movePoints(GeoMesh<polyMesh>::mesh_.points());

    return true;
}


void Foam::pointMesh::updateMesh(const mapPolyMesh& mpm)
{
    if (debug)
    {
        Pout<< "pointMesh::updateMesh(const mapPolyMesh&): "
            << "Updating for topology changes." << nl << endl;
    }
    boundary_.updateMesh();

    // Map all registered point fields
    mapFields(mpm);
}


bool Foam::pointMesh::writeObject
(
    IOstreamOption streamOpt,
    const bool writeOnProc
) const
{
    if (debug)
    {
        Pout<< "pointMesh::writeObject(IOstreamOption, const bool): "
            << "Writing to " << boundary_.objectRelPath() << endl;
    }
    return boundary_.writeObject(streamOpt, writeOnProc);
}


// ************************************************************************* //
