/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class DynamicMeshType>
Foam::simplifiedMeshes::SimplifiedDynamicFvMesh<DynamicMeshType>::
SimplifiedDynamicFvMesh
(
    const Time& runTime,
    const word& regionName
)
:
    simplifiedDynamicFvMeshBase(),
    columnFvMeshInfo(runTime, regionName),
    DynamicMeshType
    (
        IOobject
        (
            regionName,
            runTime.constant(),
            runTime,
            IOobject::NO_READ, // Do not read any existing mesh
            IOobject::NO_WRITE
        ),
        std::move(points1D_),
        std::move(faces1D_),
        std::move(owner1D_),
        std::move(neighbour1D_)
    )
{
    // Workaround to read fvSchemes and fvSolution after setting NO_READ
    // when creating the mesh
    {
        fvSchemes::readOpt(IOobject::MUST_READ);
        fvSchemes::read();
        fvSolution::readOpt(IOobject::MUST_READ);
        fvSolution::read();
    }

    // Add the patches
    addLocalPatches(*this);

    // Add the zones if constructed from mesh
    initialiseZones(*this);
}


// ************************************************************************* //
