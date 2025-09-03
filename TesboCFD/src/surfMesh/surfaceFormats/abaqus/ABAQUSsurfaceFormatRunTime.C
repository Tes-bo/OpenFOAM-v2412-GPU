/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "ABAQUSsurfaceFormat.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// Read MeshedSurface - tag as abaqus or abq
// avoid .inp (name conflict with starcd)
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    ABAQUSsurfaceFormat,
    face,
    fileExtension,
    abaqus
);

addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    ABAQUSsurfaceFormat,
    face,
    fileExtension,
    abq
);

// Write with MeshedSurfaceProxy
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    face,
    write,
    fileExtension,
    abaqus
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    face,
    write,
    fileExtension,
    abq
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    triFace,
    write,
    fileExtension,
    abaqus
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    triFace,
    write,
    fileExtension,
    abq
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    abaqus
);

addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    ABAQUSsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    abq
);

}
}

// ************************************************************************* //
