/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "NASsurfaceFormat.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// Read MeshedSurface - .bdf (Bulk Data Format) and nas (Nastran)
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    face,
    fileExtension,
    nastran
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    face,
    fileExtension,
    bdf
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    face,
    fileExtension,
    nas
);

addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    triFace,
    fileExtension,
    nastran
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    triFace,
    fileExtension,
    bdf
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    triFace,
    fileExtension,
    nas
);

addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    labelledTri,
    fileExtension,
    nastran
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    labelledTri,
    fileExtension,
    bdf
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    NASsurfaceFormat,
    labelledTri,
    fileExtension,
    nas
);


// Write MeshedSurfaceProxy
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    face,
    write,
    fileExtension,
    nastran
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    face,
    write,
    fileExtension,
    nas
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    triFace,
    write,
    fileExtension,
    nastran
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    triFace,
    write,
    fileExtension,
    nas
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    nastran
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    NASsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    nas
);

}
}

// ************************************************************************* //
