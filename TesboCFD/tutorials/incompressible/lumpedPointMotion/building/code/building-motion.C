/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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
    building-motion

Description
    Forced oscillation code for fluid-structure interface check.
    Generates position and rotation angle of each node.
    The top displacements of the target building calculated as sin() function.
    Horizontal displacements of each node interpolated to the vertical
    direction according to cantilever beam theory.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
#include "Fstream.H"
#include "unitConversion.H"
#include "foamVtkSeriesWriter.H"
#include "lumpedPointTools.H"
#include "lumpedPointState.H"
#include "lumpedPointIOMovement.H"

using namespace Foam;

//- Oscillator generator
class position_generator
{
    // Private Member Functions

        //- Calculate position/rotation at given time
        lumpedPointState calc(scalar currTime) const
        {
            // Point positions
            pointField points_(nDivisions+1, Zero);

            // Point rotations
            vectorField angles_(nDivisions+1, Zero);

            // Set node heights (z)
            forAll(points_, divi)
            {
                points_[divi].z() = (height * divi)/scalar(nDivisions);
            }

            const vector sines
            (
                Foam::sin(2*constant::mathematical::pi * currTime/period.x()),
                Foam::sin(2*constant::mathematical::pi * currTime/period.y()),
                Foam::sin(2*constant::mathematical::pi * currTime/period.z())
            );

            for (label divi = 1; divi <= nDivisions; ++divi)
            {
                const scalar zpos = points_[divi].z();
                const scalar height1 = (height - zpos);

                const scalar pos_factor =
                (
                    1.0/3.0 / pow4(height)
                  * (
                        3*pow4(height)
                      - 4*pow3(height)*(height1)
                      + pow4(height1)
                    )
                );

                const scalar ang_factor =
                (
                    1.0/3.0 / pow4(height)
                  * (
                        4*pow3(height)
                      - 4*pow3(height1)
                    )
                );

                vector here
                (
                    (amplitude.x() * sines.x() * pos_factor),
                    (amplitude.y() * sines.y() * pos_factor),
                    zpos // Z position is invariant
                );

                vector rot
                (
                    Foam::atan(amplitude.x() * sines.x() * ang_factor),
                    Foam::atan(amplitude.y() * sines.y() * ang_factor),
                    Foam::atan(amplitude.z() * sines.z() * ang_factor)
                );

                // Assign
                points_[divi] = here;

                // The x<->y swap is intentional
                angles_[divi] = vector{rot.y(), rot.x(), rot.z()};
            }

            // Return as lumpedPoint state
            return lumpedPointState{points_, angles_};
        }


public:

    // Control parameters

        // The number of oscillating nodes
        label nDivisions = 10;

        // Height of target building [m]
        scalar height = 0.5;

        // Proper period (sec)
        vector period = vector{1, 0.5, 1};

        // Amplitude
        vector amplitude = vector{0.03, 0.05, 0.3};


    // Constructors

        //- Default construct
        position_generator() = default;


    // Member Functions

        //- Calculate position/rotation at given time
        lumpedPointState state(const scalar currTime) const
        {
            return calc(currTime);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Forced oscillation code for fluid-structure interface check."
        " Generates position and rotation angle of each node."
        " The top displacements of the target building calculated as sin()"
        " function."
        " Horizontal displacements of each node interpolated to the vertical"
        " direction according to cantilever beam theory."
     );

    argList::noBanner();
    argList::noParallel();

    // Geometric
    argList::addOption
    (
        "nodes",
        "N",
        "The number of oscillating nodes (default: 10)"
    );
    argList::addOption
    (
        "height",
        "value",
        "Height of target building (m)"
    );
    argList::addOption
    (
        "period",
        "(time time time)",
        "The proper period (sec)"
    );
    argList::addOption
    (
        "amplitude",
        "(value value value)",
        "The amplitude"
    );

    // Time control
    argList::addOption
    (
        "time",
        "value",
        "The time to use"
    );
    argList::addOption
    (
        "deltaT",
        "value",
        "The time increment for multiple time loops"
    );
    argList::addOption
    (
        "nTimes",
        "value",
        "The number of time loops"
    );

    // Query, output
    argList::addBoolOption
    (
        "query",
        "Report values only and exit"
    );
    argList::addOption
    (
        "output",
        "file",
        "write to file, with header"
    );

    argList::addOption
    (
        "scale",
        "factor",
        "Scaling factor for movement (default: 1)"
    );
    argList::addOption
    (
        "visual-length",
        "len",
        "Visualization length for planes (visualized as triangles)"
    );

    // Run controls
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
    #include "setRootCase.H"


    // The oscillator
    position_generator gen;
    args.readIfPresent("nodes", gen.nDivisions);
    args.readIfPresent("height", gen.height);
    args.readIfPresent("period", gen.period);
    args.readIfPresent("amplitude", gen.amplitude);


    // Control parameters
    const bool dryrun = args.dryRun();
    const bool slave = args.found("slave");
    const bool removeLock = args.found("removeLock");

    const bool optQuery = args.found("query");
    const fileName outputFile(args.getOrDefault<fileName>("output", ""));

    const scalar relax = args.getOrDefault<scalar>("scale", 1);

    args.readIfPresent("visual-length", lumpedPointState::visLength);

    // Time parameters
    scalar currTime = args.getOrDefault<scalar>("time", 0);
    const scalar deltaT = args.getOrDefault("deltaT", 0.001);
    const label nTimes = args.getOrDefault<label>("nTimes", 1);

    // Loop handling for slave
    const bool infiniteLoop = slave && !args.found("nTimes");


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

        externalFileCoupler& coupler = movement.coupler();

        for (label timei = 0; infiniteLoop || (timei < nTimes); ++timei)
        {
            Info<< args.executable() << ": waiting for master" << endl;

            // Wait for master, but stop if status=done was seen
            if (!coupler.waitForMaster())
            {
                Info<< args.executable()
                    << ": stopping status=done was detected" << endl;
                break;
            }

            scalar timeValue = currTime;

            if (infiniteLoop)
            {
                // Get output file
                IFstream is(coupler.resolveFile(movement.outputName()));

                dictionary dict;
                is >> dict;

                timeValue = dict.get<scalar>("time");
            }

            lumpedPointState state(gen.state(timeValue));
            state.relax(relax, state0);

            // Generate input for OpenFOAM
            {
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
                    os.writeEntry("time", timeValue);
                    state.writeDict(os);
                }
            }

            Info<< args.executable()
                << ": updating state " << timei
                << " - switch to master"
                << endl;

            // Let OpenFOAM know that it can continue
            coupler.useMaster();

            currTime += deltaT;
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

        autoPtr<Time> runTimePtr;
        autoPtr<lumpedPointIOMovement> movementPtr;

        const bool oldThrowingIO = FatalIOError.throwing(true);
        const bool oldThrowingEx = FatalError.throwing(true);

        try
        {
            Info<< "Create time" << flush;

            runTimePtr = Time::New(args);

            // Create movement without a mesh
            movementPtr = lumpedPointIOMovement::New(*runTimePtr);
        }
        catch (...)
        {
            Info<< " ... failed (optional for dry-run)";
        }
        Info<< nl << endl;

        FatalError.throwing(oldThrowingEx);
        FatalIOError.throwing(oldThrowingIO);

        if (!movementPtr)
        {
            Info<< "No time, run without movement information\n" << endl;
        }

        const lumpedPointState state0(gen.state(0));

        vtk::seriesWriter stateSeries;

        for
        (
            label timei = 0, outputCount = 0;
            timei < nTimes;
            ++timei
        )
        {
            lumpedPointState state(gen.state(currTime));
            state.relax(relax, state0);

            Info<< "output [" << timei << '/' << nTimes << ']';

            // State/response = what comes back from FEM
            {
                const word outputName =
                    word::printf("state_%06d.vtp", outputCount);

                Info<< "  " << outputName;

                if (movementPtr)
                {
                    movementPtr->writeStateVTP(state, outputName);
                }
                else
                {
                    state.writeVTP(outputName);
                }
                stateSeries.append(outputCount, outputName);
            }

            Info<< endl;

            ++outputCount;
            currTime += deltaT;
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
    // Report values or generate a file
    // ----------------------------------------------------------------------

    if (optQuery || !outputFile.empty())
    {
        autoPtr<OFstream> osPtr;

        if (!outputFile.empty())
        {
            osPtr.reset(new OFstream(outputFile));
            auto& os = *osPtr;

            os.precision(8);

            // One file with everything, output using OpenFOAM syntax

            IOobject::writeBanner(os)
                << "FoamFile\n{\n"
                << "    version     " << os.version() << ";\n"
                << "    format      " << os.format() << ";\n"
                << "    class       " << "dictionary" << ";\n"
                << "    object      " << "response" << ";\n"
                << "}\n";

            IOobject::writeDivider(os) << nl;

            os << "// angles are Euler angles z-x-z (intrinsic)" << nl;
            os.writeEntry("degrees", "false");
            os << nl;
            os << "response" << nl;
            os << token::BEGIN_LIST << nl;
        }
        else
        {
            Info.stream().precision(8);
        }


        for (label timei = 0; timei < nTimes; ++timei)
        {
            lumpedPointState state(gen.state(currTime));

            if (osPtr)
            {
                // Report position/angle
                auto& os = *osPtr;

                os.beginBlock();

                os.writeEntry("time", currTime);

                state.writeDict(os);

                os.endBlock();
            }
            else
            {
                // Report position/angle
                auto& os = Info.stream();

                os.writeEntry("time", currTime);
                state.writeDict(os);
            }

            currTime += deltaT;
        }


        if (osPtr)
        {
            auto& os = *osPtr;

            os << token::END_LIST;
            os.endEntry();

            IOobject::writeEndDivider(os);

            Info<< "\nEnd\n" << endl;
        }

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
        timei < nTimes;
        ++timei
    )
    {
        lumpedPointState state(gen.state(currTime));

        state += movement.origin();
        movement.scalePoints(state);
        state.relax(relax, state0);

        Info<< "output [" << timei << '/' << nTimes << ']';

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
        currTime += deltaT;
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
