/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

#include "X3DsurfaceFormatCore.H"
#include "Ostream.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::X3DsurfaceFormatCore::writeHeader
(
    Ostream& os
)
{
    os  <<
        "<?xml version='1.0' encoding='UTF-8'?>\n"
        "<!DOCTYPE X3D PUBLIC \"ISO//Web3D//DTD X3D 3.0//EN\" "
        "\"http://www.web3d.org/specifications/x3d-3.0.dtd\">\n"
        "<X3D\n"
        "  version='3.0'\n"
        "  profile='Immersive'\n"
        "  xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance'\n"
        "  xsd:noNamespaceSchemaLocation="
        "'http://www.web3d.org/specifications/x3d-3.0.xsd'\n"
        "  >\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::writeFooter
(
    Ostream& os
)
{
    os  <<
        "</X3D>\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::beginGroup
(
    Ostream& os
)
{
    os  <<
        "<Group>\n"
        " <Shape>\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::endGroup
(
    Ostream& os
)
{
    os  <<
        "  </Shape>\n"
        " </Group>\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::writeAppearance
(
    Ostream& os
)
{
    os  <<
        "  <Appearance>\n"
        "   <Material"
        " ambientIntensity='0'"
        " diffuseColor='1 1 1'" // Default: '0.8 0.8 0.8'
        // Default: " emissiveColor='0 0 0'"
        // Default: " specularColor='0 0 0'"
        " shininess='0.8'"      // Default: 0.2
        " transparency='0'"
        " />\n"           // Material
        "  </Appearance>\n";
}


void Foam::fileFormats::X3DsurfaceFormatCore::writePoints
(
    Ostream& os,
    const UList<point>& pts
)
{
    os  <<
        "    <Coordinate point='\n";

    for (const point& p : pts)
    {
        os  << p.x() << ' ' << p.y() << ' ' << p.z() << ',' << nl;
    }

    os  <<
        "' />\n";
}


// ************************************************************************* //
