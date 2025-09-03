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
    surfaceOrient

Group
    grpSurfaceUtilities

Description
    Set normal consistent with respect to a user provided 'outside' point.
    If the -inside option is used the point is considered inside.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurfaceSearch.H"
#include "orientedSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Set face normals consistent with a user-provided 'outside' point"
    );

    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("point", "The visible 'outside' point");
    argList::addArgument("output", "The output surface file");

    argList::addBoolOption
    (
        "inside",
        "Treat provided point as being inside"
    );
    argList::addBoolOption
    (
        "usePierceTest",
        "Determine orientation by counting number of intersections"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Input geometry scaling factor"
    );

    argList args(argc, argv);

    const auto surfFileName = args.get<fileName>(1);
    const auto visiblePoint = args.get<point>(2);
    const auto outFileName  = args.get<fileName>(3);

    const bool orientInside = args.found("inside");
    const bool usePierceTest = args.found("usePierceTest");

    Info<< "Reading surface from " << surfFileName << nl
        << "Orienting surface such that visiblePoint " << visiblePoint
        << " is ";

    if (orientInside)
    {
        Info<< "inside" << endl;
    }
    else
    {
        Info<< "outside" << endl;
    }

    const scalar scaling = args.getOrDefault<scalar>("scale", -1);
    if (scaling > 0)
    {
        Info<< "Input scaling: " << scaling << nl;
    }

    // Load surface
    triSurface surf(surfFileName, scaling);

    bool anyFlipped = false;

    if (usePierceTest)
    {
        triSurfaceSearch surfSearches(surf);

        anyFlipped = orientedSurface::orient
        (
            surf,
            surfSearches,
            visiblePoint,
           !orientInside
        );
    }
    else
    {
        anyFlipped = orientedSurface::orient
        (
            surf,
            visiblePoint,
           !orientInside
        );
    }

    if (anyFlipped)
    {
        Info<< "Flipped orientation of (part of) surface." << nl;
    }
    else
    {
        Info<< "Did not flip orientation of any triangle of surface." << nl;
    }

    Info<< "Writing new surface to " << outFileName << endl;

    surf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
