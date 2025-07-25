const dictionary& pimple = mesh.solutionDict().subDict("PIMPLE");

int nNonOrthCorr =
    pimple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
