Info<< "Total cylinder mass: " << fvc::domainIntegrate(rho).value() << endl;

autoPtr<OFstream> logSummaryFile;

if (Pstream::master())
{
    logSummaryFile.reset
    (
        new OFstream
        (
            runTime.globalPath()
           /("logSummary." + runTime.timeName() + ".dat")
        )
    );

    logSummaryFile()
        << "# CA" << "      p" << "        T" << "       u'" << endl;
}
