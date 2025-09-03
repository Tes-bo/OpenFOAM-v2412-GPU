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

#include "STARCDsurfaceFormat.H"
#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// Read MeshedSurface
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STARCDsurfaceFormat,
    face,
    fileExtension,
    starcd
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STARCDsurfaceFormat,
    triFace,
    fileExtension,
    starcd
);
addNamedTemplatedToRunTimeSelectionTable
(
    MeshedSurface,
    STARCDsurfaceFormat,
    labelledTri,
    fileExtension,
    starcd
);

// Write MeshedSurfaceProxy
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STARCDsurfaceFormat,
    face,
    write,
    fileExtension,
    starcd
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STARCDsurfaceFormat,
    triFace,
    write,
    fileExtension,
    starcd
);
addNamedTemplatedToMemberFunctionSelectionTable
(
    MeshedSurfaceProxy,
    STARCDsurfaceFormat,
    labelledTri,
    write,
    fileExtension,
    starcd
);

}
}

// ************************************************************************* //
