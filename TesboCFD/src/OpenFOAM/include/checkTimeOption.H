// Deprecated include (2024-05) - prefer timeSelector
// ----------------------------------------------------------------------------

// Check -time and -latestTime options
if (args.found("latestTime"))
{
    startTime = Times.size() - 1;
}
else if (args.found("time"))
{
    Foam::scalar timeValue = args.get<Foam::scalar>("time");

    startTime = Foam::Time::findClosestTimeIndex(Times, timeValue);
}


// ************************************************************************* //
