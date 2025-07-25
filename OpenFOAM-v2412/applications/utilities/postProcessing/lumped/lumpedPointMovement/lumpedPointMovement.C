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
    lumpedPointMovement

Description
    This utility can be used to produce VTK files to visualize the response
    points/rotations and the corresponding movement of the building surfaces.

    Uses the tabulated responses from the specified file.
    Optionally, it can also be used to a dummy responder for the
    externalFileCoupler logic, which makes it useful as a debugging facility
    as well demonstrating how an external application could communicate
    with the lumpedPointMovement point-patch boundary condition.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "OFstream.H"
#include "foamVtkSeriesWriter.H"
#include "lumpedPointTools.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;


inline List<lumpedPointStateTuple> getResponseTable
(
    const fileName& file,
    const lumpedPointState& state0
)
{
    return lumpedPointTools::lumpedPointStates
    (
        file,
        state0.rotationOrder(),
        state0.degrees()
    );
}


void echoTableLimits
(
    const List<lumpedPointStateTuple>& tbl,
    const label span,
    const label maxOut
)
{
    Info<< "Using response table with " << tbl.size() << " entries" << nl;

    if (span)
    {
        Info<< "Increment input by " << span << nl;
    }

    if (maxOut)
    {
        Info<< "Stopping after " << maxOut << " outputs" << nl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Visualize lumpedPoint movements or provide a slave responder"
        " for diagnostic purposes."
    );

    argList::noFunctionObjects();  // Never use function objects
    argList::addOption
    (
        "max",
        "N",
        "Maximum number of outputs"
    );
    argList::addOption
    (
        "span",
        "N",
        "Increment each input by N (default: 1)"
    );
    argList::addOption
    (
        "scale",
        "factor",
        "Relaxation/scaling factor for movement (default: 1)"
    );
    argList::addOption
    (
        "visual-length",
        "len",
        "Visualization length for planes (visualized as triangles)"
    );
    argList::addDryRunOption
    (
        "Test movement without a mesh"
    );
    argList::addBoolOption
    (
        "removeLock",
        "Remove lock-file on termination of slave"
    );
    argList::addBoolOption
    (
        "slave",
        "Invoke as a slave responder for testing"
    );
    argList::addArgument("responseFile");

    #include "setRootCase.H"

    const label maxOut = Foam::max(0, args.getOrDefault<label>("max", 0));
    const label span   = Foam::max(1, args.getOrDefault<label>("span", 1));

    // Control parameters
    const bool dryrun = args.dryRun();
    const bool slave = args.found("slave");
    const bool removeLock = args.found("removeLock");

    const scalar relax = args.getOrDefault<scalar>("scale", 1);

    args.readIfPresent("visual-length", lumpedPointState::visLength);

    const auto responseFile = args.get<fileName>(1);

    // ----------------------------------------------------------------------
    // Slave mode
    // ----------------------------------------------------------------------

    if (slave)
    {
        Info<< "Running as slave responder" << endl;

        if (Pstream::parRun())
        {
            FatalErrorInFunction
                << "Running as slave responder is not permitted in parallel"
                << nl
                << exit(FatalError);
        }

        #include "createTime.H"

        // Create movement without a mesh
        autoPtr<lumpedPointIOMovement> movementPtr =
            lumpedPointIOMovement::New(runTime);

        if (!movementPtr)
        {
            Info<< "No valid movement found" << endl;
            return 1;
        }
        auto& movement = *movementPtr;

        // Reference state0
        const lumpedPointState& state0 = movement.state0();

        List<lumpedPointStateTuple> responseTable =
            getResponseTable(responseFile, state0);

        echoTableLimits(responseTable, span, maxOut);

        if (dryrun)
        {
            Info<< "dry-run: response table with " << responseTable.size()
                << " entries" << nl
                << "\nEnd\n" << endl;
            return 0;
        }

        externalFileCoupler& coupler = movement.coupler();

        for
        (
            label timei = 0, outputCount = 0;
            timei < responseTable.size();
            timei += span
        )
        {
            Info<< args.executable() << ": waiting for master" << endl;

            // Wait for master, but stop if status=done was seen
            if (!coupler.waitForMaster())
            {
                Info<< args.executable()
                    << ": stopping status=done was detected" << endl;
                break;
            }

            lumpedPointState state = responseTable[timei].second();
            state.relax(relax, state0);

            // Generate input for OpenFOAM
            OFstream os(coupler.resolveFile(movement.inputName()));
            if
            (
                movement.inputFormat()
             == lumpedPointState::inputFormatType::PLAIN
            )
            {
                state.writePlain(os);
            }
            else
            {
                os.writeEntry("time", responseTable[timei].first());
                state.writeDict(os);
            }

            Info<< args.executable()
                << ": updated to state " << timei
                << " - switch to master"
                << endl;

            // Let OpenFOAM know that it can continue
            coupler.useMaster();

            ++outputCount;

            if (maxOut && outputCount >= maxOut)
            {
                Info<< args.executable()
                    << ": stopping after " << maxOut << " outputs" << endl;
                break;
            }
        }

        if (removeLock)
        {
            Info<< args.executable() << ": removing lock file" << endl;
            coupler.useSlave();  // This removes the lock-file
        }

        Info<< args.executable() << ": finishing" << nl;

        Info<< "\nEnd\n" << endl;
        return 0;
    }


    // ----------------------------------------------------------------------
    // dry-run
    // ----------------------------------------------------------------------

    if (dryrun)
    {
        Info<< "dry-run: creating states only" << nl;

        #include "createTime.H"

        // Create movement without a mesh
        autoPtr<lumpedPointIOMovement> movementPtr =
            lumpedPointIOMovement::New(runTime);

        if (!movementPtr)
        {
            Info<< "No valid movement found" << endl;
            return 1;
        }
        auto& movement = *movementPtr;

        // Reference state0
        const lumpedPointState& state0 = movement.state0();

        List<lumpedPointStateTuple> responseTable =
            getResponseTable(responseFile, state0);

        echoTableLimits(responseTable, span, maxOut);


        vtk::seriesWriter stateSeries;

        for
        (
            label timei = 0, outputCount = 0;
            timei < responseTable.size();
            timei += span
        )
        {
            lumpedPointState state = responseTable[timei].second();

            state += movement.origin();
            movement.scalePoints(state);
            state.relax(relax, state0);

            Info<< "output [" << timei << '/' << responseTable.size() << ']';

            // State/response = what comes back from FEM
            {
                const word outputName =
                    word::printf("state_%06d.vtp", outputCount);

                Info<< "  " << outputName;

                movement.writeStateVTP(state, outputName);
                stateSeries.append(outputCount, outputName);
            }

            Info<< endl;

            ++outputCount;

            if (maxOut && outputCount >= maxOut)
            {
                Info<< "Max output " << maxOut << " ... stopping" << endl;
                break;
            }
        }

        // Write file series

        if (stateSeries.size())
        {
            Info<< nl << "write state.vtp.series" << nl;
            stateSeries.write("state.vtp");
        }

        Info<< "\nEnd\n" << endl;
        return 0;
    }


    // ----------------------------------------------------------------------
    // test patch movement
    // ----------------------------------------------------------------------

    #include "createTime.H"

    runTime.setTime(instant(0, runTime.constant()), 0);

    #include "createNamedMesh.H"

    // Create movement with mesh
    autoPtr<lumpedPointIOMovement> movementPtr =
        lumpedPointIOMovement::New(mesh);

    if (!movementPtr)
    {
        Info<< "No valid movement found" << endl;
        return 1;
    }
    auto& movement = *movementPtr;

    // Reference state0
    const lumpedPointState& state0 = movement.state0();

    List<lumpedPointStateTuple> responseTable =
        getResponseTable(responseFile, state0);

    echoTableLimits(responseTable, span, maxOut);

    pointIOField points0(lumpedPointTools::points0Field(mesh));

    const label nPatches = lumpedPointTools::setPatchControls(mesh, points0);
    if (!nPatches)
    {
        Info<< "No point patches with lumped movement found" << endl;
        return 2;
    }

    Info<< "Lumped point patch controls set on "
        << nPatches << " patches" << nl;

    lumpedPointTools::setInterpolators(mesh, points0);


    // Output vtk file series
    vtk::seriesWriter stateSeries;
    vtk::seriesWriter geomSeries;

    // Initial geometry
    movement.writeVTP("geom_init.vtp", state0, mesh, points0);

    lumpedPointTools::setInterpolators(mesh);

    for
    (
        label timei = 0, outputCount = 0;
        timei < responseTable.size();
        timei += span
    )
    {
        lumpedPointState state = responseTable[timei].second();

        state += movement.origin();
        movement.scalePoints(state);
        state.relax(relax, state0);

        Info<< "output [" << timei << '/' << responseTable.size() << ']';

        // State/response = what comes back from FEM
        {
            const word outputName =
                word::printf("state_%06d.vtp", outputCount);

            Info<< "  " << outputName;

            movement.writeStateVTP(state, outputName);
            stateSeries.append(outputCount, outputName);
        }

        {
            const word outputName =
                word::printf("geom_%06d.vtp", outputCount);

            Info<< "  " << outputName;

            movement.writeVTP(outputName, state, mesh, points0);
            geomSeries.append(outputCount, outputName);
        }

        Info<< endl;

        ++outputCount;

        if (maxOut && outputCount >= maxOut)
        {
            Info<< "Max output " << maxOut << " ... stopping" << endl;
            break;
        }
    }


    // Write file series

    if (geomSeries.size())
    {
        Info<< nl << "write geom.vtp.series" << nl;
        geomSeries.write("geom.vtp");
    }
    if (stateSeries.size())
    {
        Info<< nl << "write state.vtp.series" << nl;
        stateSeries.write("state.vtp");
    }

    Info<< "\nEnd\n" << endl;

    return 0;
}


// ************************************************************************* //
