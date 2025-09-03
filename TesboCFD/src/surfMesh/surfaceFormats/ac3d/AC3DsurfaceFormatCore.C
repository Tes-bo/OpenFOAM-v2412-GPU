/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "AC3DsurfaceFormatCore.H"
#include "clock.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// Define 8 standard colours as r,g,b components
static float colourMap[] =
{
    1, 1, 1,
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 0,
    0, 1, 1,
    1, 0, 1,
    0.5, 0.5, 1
};


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

bool Foam::fileFormats::AC3DsurfaceFormatCore::readCmd
(
    IFstream& is,
    string& cmd,
    string& args
)
{
    if (is.good())
    {
        string line;
        is.getLine(line);

        const auto space = line.find(' ');

        if (space && space != string::npos)
        {
            cmd  = line.substr(0, space);
            args = line.substr(space+1);

            return true;
        }
    }

    return false;
}


bool Foam::fileFormats::AC3DsurfaceFormatCore::cueTo
(
    IFstream& is,
    const string& cmd,
    string& args
)
{
    string line;
    while (is.good())
    {
        is.getLine(line);

        const auto space = line.find(' ');

        if (space && space != string::npos && cmd == line.substr(0, space))
        {
            args = line.substr(space+1);

            return true;
        }
    }

    return false;
}


Foam::string Foam::fileFormats::AC3DsurfaceFormatCore::cueToOrDie
(
    IFstream& is,
    const string& cmd,
    const string& errorMsg
)
{
    string args;
    if (!cueTo(is, cmd, args))
    {
        FatalErrorInFunction
            << "Cannot find command " << cmd
            << " " << errorMsg
            << exit(FatalError);
    }

    return args;
}


void Foam::fileFormats::AC3DsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const UList<surfZone>& zones
)
{
    // Write with zones as separate objects under "world" object.
    // Header is taken over from sample file.
    // Defines separate materials for all zones. Recycle colours.

    // Write header. Define materials.
    os  << "AC3Db" << nl;

    forAll(zones, zonei)
    {
        const label colourI = zonei % 8;
        const label colourCompI = 3 * colourI;

        os  << "MATERIAL \"" << zones[zonei].name() << "Mat\" rgb "
            << colourMap[colourCompI] << ' ' << colourMap[colourCompI+1]
            << ' ' << colourMap[colourCompI+2]
            << "  amb 0.2 0.2 0.2  emis 0 0 0  spec 0.5 0.5 0.5  shi 10"
            << "  trans 0"
            << nl;
    }

    os  << "OBJECT world" << nl
        << "kids " << zones.size() << endl;
}


// ************************************************************************* //
