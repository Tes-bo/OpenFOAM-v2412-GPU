{
    const scalar meanP = p.weightedAverage(mesh.V()).value();
    const scalar meanT = thermo.T().weightedAverage(mesh.V()).value();
    const scalar meanUp =
        (sqrt((2.0/3.0)*turbulence->k()))().weightedAverage(mesh.V()).value();
    const scalar meanB = b.weightedAverage(mesh.V()).value();

    Info<< "Mean pressure:" << meanP << nl
        << "Mean temperature:" << meanT << nl
        << "Mean u':" << meanUp << endl;

    if (Pstream::master())
    {
        logSummaryFile()
            << runTime.theta() << tab
            << meanP << tab
            << meanT << tab
            << meanUp << tab
            << 1 - meanB
            << endl;
    }
}
