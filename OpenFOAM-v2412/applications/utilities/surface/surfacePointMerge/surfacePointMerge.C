/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    surfacePointMerge

Group
    grpSurfaceUtilities

Description
    Merges points on surface if they are within absolute distance.
    Since absolute distance use with care!

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Merge points on surface if they are within absolute distance [m]."
    );
    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("distance", "The merge distance");
    argList::addArgument("output", "The output surface file");

    argList::addOption
    (
        "scale",
        "factor",
        "Input geometry scaling factor"
    );

    argList args(argc, argv);

    const auto surfFileName = args.get<fileName>(1);
    const auto mergeTol = args.get<scalar>(2);
    const auto outFileName = args.get<fileName>(3);

    const scalar scaling = args.getOrDefault<scalar>("scale", -1);

    Info<< "Reading surface from " << surfFileName << " ..." << nl
        << "Merging points within " << mergeTol << " metre." << nl;
    if (scaling > 0)
    {
        Info<< "input scaling " << scaling << nl;
    }

    const triSurface surf1(surfFileName, scaling);

    Info<< "Original surface:" << nl;
    surf1.writeStats(Info);

    triSurface cleanSurf(surf1);

    while (true)
    {
        const label nOldVert = cleanSurf.nPoints();

        cleanSurf = triSurfaceTools::mergePoints(cleanSurf, mergeTol);

        Info<< "After merging points:" << endl;

        cleanSurf.writeStats(Info);

        if (nOldVert == cleanSurf.nPoints())
        {
            break;
        }
    }

    cleanSurf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
