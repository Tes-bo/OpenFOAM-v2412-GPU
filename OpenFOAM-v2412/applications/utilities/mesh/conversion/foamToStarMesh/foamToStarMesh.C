/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Application
    foamToStarMesh

Group
    grpMeshConversionUtilities

Description
    Write an OpenFOAM mesh in STARCD/PROSTAR (v4) bnd/cel/vrt format.

Usage
    \b foamToStarMesh [OPTION]

    Options:
      - \par -noBnd
        Suppress writing the \c .bnd file

      - \par -scale \<factor\>
        Specify an alternative geometry scaling factor.
        The default is \b 1000 (scale \em [m] to \em [mm]).

Note
    The cellTable information available in the files
    \c constant/cellTable and \c constant/polyMesh/cellTableId
    will be used if available. Otherwise the cellZones are used when
    creating the cellTable information.

See also
    Foam::cellTable, Foam::meshWriter and Foam::fileFormats::STARCDMeshWriter

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "polyMesh.H"
#include "STARCDMeshWriter.H"
#include "IOdictionary.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Write an OpenFOAM mesh in STARCD/PROSTAR (v4) bnd/cel/vrt format"
    );
    argList::noParallel();
    timeSelector::addOptions();

    argList::addOption
    (
        "scale",
        "factor",
        "Geometry scaling factor - default is 1000 ([m] to [mm])"
    );
    argList::addBoolOption
    (
        "noBnd",
        "Suppress writing a boundary (.bnd) file"
    );

    #include "setRootCase.H"
    #include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

    fileName exportName = meshWriter::defaultMeshName;
    if (args.found("case"))
    {
        exportName += '-' + args.globalCaseName();
    }

    // Default rescale from [m] to [mm]
    const scalar scaleFactor = args.getOrDefault<scalar>("scale", 1000);
    const bool  writeBndFile = !args.found("noBnd");

    #include "createPolyMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        #include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            fileFormats::STARCDMeshWriter writer
            (
                mesh,
                scaleFactor,
                writeBndFile
            );

            fileName meshName(exportName);
            if (state != polyMesh::UNCHANGED)
            {
                meshName += '_' + runTime.timeName();
            }

            writer.write(meshName);
        }

        Info<< nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
