#include "createControl.H"

int nCorr = stressControl.getOrDefault<int>("nCorrectors", 1);

scalar convergenceTolerance(stressControl.get<scalar>("D"));
