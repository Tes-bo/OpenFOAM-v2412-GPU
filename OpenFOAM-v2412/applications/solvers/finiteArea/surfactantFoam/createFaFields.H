Info<< "Reading field Cs" << endl;
areaScalarField Cs
(
    IOobject
    (
        "Cs",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    aMesh
);

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


Info<< "Reading diffusivity D\n" << endl;

dimensionedScalar Ds("Ds", dimViscosity, transportProperties);

areaVectorField Us
(
    IOobject
    (
        "Us",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    aMesh
);


edgeScalarField phis
(
    IOobject
    (
        "phis",
        runTime.timeName(),
        aMesh.thisDb(),
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearEdgeInterpolate(Us) & aMesh.Le()
);
