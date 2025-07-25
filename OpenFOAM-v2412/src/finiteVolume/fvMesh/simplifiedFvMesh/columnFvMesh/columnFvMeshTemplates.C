/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018 OpenCFD Ltd.
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
#include "polyMesh.H"

template<class ZoneMeshType>
void Foam::simplifiedMeshes::columnFvMeshInfo::initialiseZone
(
    const word& zoneTypeName,
    const fileName& instance,
    ZoneMeshType& zoneType
)
{
    const wordList zoneNames
    (
        ZoneMeshType
        (
            IOobject
            (
                zoneTypeName + "Zones",
                instance,
                polyMesh::meshSubDir,
                zoneType.mesh(),
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            zoneType.mesh()
        ).names()
    );

    ZoneMeshType::disallowGenericZones = 1;
    for (const word& name : zoneNames)
    {
        // Insert empty zone
        (void)zoneType[name];
    }
    ZoneMeshType::disallowGenericZones = 0;
}


// ************************************************************************* //
