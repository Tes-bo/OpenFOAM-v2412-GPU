const dictionary& stressControl = mesh.solutionDict().subDict("stressAnalysis");

bool compactNormalStress(stressControl.get<bool>("compactNormalStress"));
