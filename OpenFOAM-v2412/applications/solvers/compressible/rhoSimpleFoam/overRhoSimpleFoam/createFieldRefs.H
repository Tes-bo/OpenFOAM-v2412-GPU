const volScalarField& psi = thermo.psi();

bool adjustFringe
(
    simple.dict().getOrDefault("oversetAdjustPhi", false)
);
