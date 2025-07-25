#include "readDyMControls.H"

correctPhi = pimple.dict().getOrDefault("correctPhi", true);
maxAcousticCo = runTime.controlDict().get<scalar>("maxAcousticCo");
