    const wordList fluidNames(rp["fluid"]);

    PtrList<fvMesh> fluidRegions(fluidNames.size());

    forAll(fluidNames, i)
    {
        if (args.dryRun() || args.found("dry-run-write"))
        {
            Info
                << "Operating in 'dry-run' mode:"
                << " case will run for 1 time step.  "
                << "All checks assumed OK on a clean exit" << endl;

            FieldBase::allowConstructFromLargerSize = true;

            // Create a simplified 1D mesh and attempt to re-create boundary
            // conditions
            fluidRegions.set
            (
                i,
                new simplifiedMeshes::columnFvMesh(runTime, fluidNames[i])
            );

            // Stop after 1 iteration of the simplified mesh

            if (args.found("dry-run-write"))
            {
                // Using saWriteNow triggers function objects execute(), write()
                runTime.stopAt(Time::saWriteNow);

                // Make sure mesh gets output to the current time (since
                // instance no longer constant)
                fluidRegions[i].setInstance(runTime.timeName());
            }
            else
            {
                // Using saNoWriteNow triggers function objects execute(),
                // but not write()
                runTime.stopAt(Time::saNoWriteNow);
            }

            functionObject::outputPrefix = "postProcessing-dry-run";
        }
        else
        {
            Info<< "Create fluid mesh for region " << fluidNames[i]
                << " for time = " << runTime.timeName() << nl << endl;

            fluidRegions.set
            (
                i,
                new fvMesh
                (
                    IOobject
                    (
                        fluidNames[i],
                        runTime.timeName(),
                        runTime,
                        IOobject::MUST_READ
                    )
                )
            );
        }
    }
