// Volume-to surface mapping object
const volSurfaceMapping vsm(aMesh);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector(dimVelocity, Zero)
);


volScalarField H
(
    IOobject
    (
        "H",
        runTime.timeName(),
        mesh.thisDb(),
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimLength, Zero)
);
