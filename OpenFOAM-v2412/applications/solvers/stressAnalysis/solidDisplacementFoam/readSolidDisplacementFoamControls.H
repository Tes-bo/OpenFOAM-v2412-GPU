nCorr = stressControl.getOrDefault<int>("nCorrectors", 1);
convergenceTolerance = stressControl.get<scalar>("D");
compactNormalStress = stressControl.get<bool>("compactNormalStress");
