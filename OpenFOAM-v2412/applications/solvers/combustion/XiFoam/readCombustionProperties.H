    Info<< "Reading combustion properties\n" << endl;

    IOdictionary combustionProperties
    (
        IOobject
        (
            "combustionProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    const word SuModel
    (
        combustionProperties.get<word>("SuModel")
    );

    dimensionedScalar sigmaExt
    (
        "sigmaExt", dimless/dimTime, combustionProperties
    );

    const word XiModel
    (
        combustionProperties.get<word>("XiModel")
    );

    dimensionedScalar XiCoef
    (
        "XiCoef", dimless, combustionProperties
    );

    dimensionedScalar XiShapeCoef
    (
        "XiShapeCoef", dimless, combustionProperties
    );

    dimensionedScalar uPrimeCoef
    (
        "uPrimeCoef", dimless, combustionProperties
    );

    ignition ign(combustionProperties, runTime, mesh);
