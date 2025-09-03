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
    surfaceRefineRedGreen

Group
    grpSurfaceUtilities

Description
    Refine by splitting all three edges of triangle ('red' refinement).

    Neighbouring triangles (which are not marked for refinement get split
    in half ('green' refinement).

    Reference:
    \verbatim
    R. Verfuerth, "A review of a posteriori
    error estimation and adaptive mesh refinement techniques",
    Wiley-Teubner, 1996)
    \endverbatim

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "triSurfaceTools.H"
#include "argList.H"
#include "OFstream.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Refine by splitting all three edges of triangle"
    );
    argList::noParallel();
    argList::addArgument("input", "The input surface file");
    argList::addArgument("output", "The output surface file");
    argList::addOption
    (
        "steps",
        "N",
        "Number of refinement steps (default: 1)"
    );
    argList args(argc, argv);

    const auto surfFileName = args.get<fileName>(1);
    const auto outFileName = args.get<fileName>(2);

    Info<< "Reading surface from " << surfFileName << " ..." << endl;

    triSurface surf(surfFileName);

    Info<< "Original surface:" << nl
        << "    triangles     :" << surf.size() << nl
        << "    vertices(used):" << surf.nPoints() << endl;


    const label nsteps =
        args.getCheckOrDefault<label>("steps", 1, labelMinMax::ge(1));


    Info<< "Refining " << nsteps << " times" << flush;

    for (label step = 0; step < nsteps; ++step)
    {
        Info<< '.' << flush;

        surf = triSurfaceTools::redGreenRefine
        (
            surf,
            identity(surf.size())  // Refine all
        );
    }
    Info<< endl;

    Info<< "Refined surface:" << nl
        << "    triangles     :" << surf.size() << nl
        << "    vertices(used):" << surf.nPoints() << endl;

    Info<< nl
        << "Writing refined surface to " << outFileName << " ..." << endl;

    surf.write(outFileName);

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
