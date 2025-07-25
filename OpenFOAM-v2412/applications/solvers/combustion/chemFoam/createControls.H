bool adjustTimeStep(runTime.controlDict().get<bool>("adjustTimeStep"));

scalar maxDeltaT(runTime.controlDict().get<scalar>("maxDeltaT"));
