const word kinematicCloudName
(
    args.getOrDefault<word>("cloud", "kinematicCloud")
);

Info<< "Constructing kinematicCloud " << kinematicCloudName << endl;

basicKinematicCloud parcels
(
    kinematicCloudName,
    rhoInf,
    U,
    muc,
    g
);

