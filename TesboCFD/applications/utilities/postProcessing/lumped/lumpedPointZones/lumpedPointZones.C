/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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
    lumpedPointZones

Description
    Produce a VTK PolyData file \c lumpedPointZones.vtp in which the
    segmentation of the pressure integration zones can be visualized
    for diagnostic purposes. Does not use external coupling.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "timeSelector.H"

#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"
#include "fvMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Create lumpedPointZones.vtp to verify the segmentation of"
        " pressure integration zones used by lumpedPoint BC."
    );

    argList::noFunctionObjects();  // Never use function objects

    argList::addDryRunOption
    (
        "Test initial lumped points state without a mesh"
    );
    argList::addOption
    (
        "visual-length",
        "len",
        "Visualization length for planes (visualized as triangles)"
    );

    argList::addBoolOption
    (
        "no-interpolate",
        "Suppress calculation/display of point interpolators"
    );

    argList::addVerboseOption();

    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"

    const bool noInterpolate = args.found("no-interpolate");

    args.readIfPresent("visual-length", lumpedPointState::visLength);

    if (args.dryRun())
    {
        // Create without a mesh
        autoPtr<lumpedPointIOMovement> movement =
            lumpedPointIOMovement::New(runTime);

        if (!movement)
        {
            Info<< "No valid movement found" << endl;
            return 1;
        }

        const word outputName("state.vtp");

        Info<< "dry-run: writing " << outputName << nl;

        movement().writeStateVTP(movement().state0(), outputName);

        Info<< "\nEnd\n" << endl;

        return 0;
    }


    runTime.setTime(instant(0, runTime.constant()), 0);

    #include "createNamedMesh.H"

    autoPtr<lumpedPointIOMovement> movement = lumpedPointIOMovement::New(mesh);

    if (!movement)
    {
        Info<< "No valid movement found" << endl;
        return 1;
    }

    // Initial positions/rotations
    movement().writeStateVTP("state.vtp");

    pointIOField points0(lumpedPointTools::points0Field(mesh));

    const label nPatches = lumpedPointTools::setPatchControls(mesh, points0);
    if (!nPatches)
    {
        Info<< "No point patches with lumped movement found" << endl;
        return 2;
    }

    Info<<"Lumped point patch controls set on "
        << nPatches << " patches" << nl;

    Info<<"Areas per point: " << flatOutput(movement().areas(mesh)) << nl;

    if (noInterpolate)
    {
        // Initial geometry, with zones
        movement().writeZonesVTP("lumpedPointZones.vtp", mesh, points0);
    }
    else
    {
        lumpedPointTools::setInterpolators(mesh, points0);

        // Initial geometry, with zones and interpolations
        movement().writeVTP("lumpedPointZones.vtp", mesh, points0);
    }

    Info<< nl
        << "wrote 'state.vtp' (reference state)" << nl
        << "wrote 'lumpedPointZones.vtp'" << nl
        << "\nEnd\n" << endl;

    return 0;
}

// ************************************************************************* //
