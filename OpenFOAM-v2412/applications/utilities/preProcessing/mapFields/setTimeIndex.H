    instantList sourceTimes = runTimeSource.times();
    label sourceTimeIndex = runTimeSource.timeIndex();
    if (args.found("sourceTime"))
    {
        if (args["sourceTime"] == "latestTime")
        {
            sourceTimeIndex = sourceTimes.size() - 1;
        }
        else
        {
            sourceTimeIndex = Time::findClosestTimeIndex
            (
                sourceTimes,
                args.get<scalar>("sourceTime")
            );
        }
    }
    else
    {
        sourceTimeIndex = Time::findClosestTimeIndex
        (
            sourceTimes,
            runTimeTarget.time().value()
        );
    }

    runTimeSource.setTime(sourceTimes[sourceTimeIndex], sourceTimeIndex);

    Info<< "\nSource time: " << runTimeSource.value()
        << "\nTarget time: " << runTimeTarget.value()
        << endl;
