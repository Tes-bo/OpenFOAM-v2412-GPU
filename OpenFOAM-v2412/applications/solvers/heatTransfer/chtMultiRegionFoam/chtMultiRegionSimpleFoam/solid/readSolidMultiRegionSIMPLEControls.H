const dictionary& simple = mesh.solutionDict().subDict("SIMPLE");

const int nNonOrthCorr =
    simple.getOrDefault<int>("nNonOrthogonalCorrectors", 0);
