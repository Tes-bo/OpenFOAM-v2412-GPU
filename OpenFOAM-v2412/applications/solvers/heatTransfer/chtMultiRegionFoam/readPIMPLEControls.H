    // We do not have a top-level mesh. Construct the fvSolution for
    // the runTime instead.
    fvSolution solutionDict(runTime);

    const dictionary& pimple = solutionDict.subDict("PIMPLE");

    const int nOuterCorr =
        pimple.getOrDefault<int>("nOuterCorrectors", 1);
