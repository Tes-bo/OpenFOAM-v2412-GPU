const dictionary& potentialFlow
(
    mesh.solutionDict().subDict("potentialFlow")
);

const int nNonOrthCorr
(
    potentialFlow.getOrDefault<int>("nNonOrthogonalCorrectors", 0)
);
