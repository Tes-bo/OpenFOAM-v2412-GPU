IOdictionary additionalControlsDict
(
    IOobject
    (
        "additionalControls",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

bool solvePrimaryRegion
(
    additionalControlsDict.getOrDefault("solvePrimaryRegion", true)
);

bool solvePyrolysisRegion
(
    additionalControlsDict.getOrDefault("solvePyrolysisRegion", true)
);

scalar maxDi = pyrolysis.maxDiff();
