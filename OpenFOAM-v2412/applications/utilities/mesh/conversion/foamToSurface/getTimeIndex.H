// Read time index from */uniform/time, but treat 0 and constant specially

    word timeName = "0";

    if
    (
        runTime.timeName() != runTime.constant()
     && runTime.timeName() != "0"
    )
    {
        IOobject io
        (
            "time",
            runTime.timeName(),
            "uniform",
            runTime,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        );

        if (io.typeHeaderOk<IOdictionary>(true))
        {
            const label timeIndex = IOdictionary(io).get<label>("index");
            timeName = Foam::name(timeIndex);
        }
        else
        {
            timeName = runTime.timeName();
            // Info<< "skip ... missing entry " << io.objectPath() << endl;
            // continue;
        }
    }

    Info<< "\nTime [" << timeName << "] = " << runTime.timeName() << nl;
