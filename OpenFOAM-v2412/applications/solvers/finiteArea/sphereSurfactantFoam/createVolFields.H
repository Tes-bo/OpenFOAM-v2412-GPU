    // Volume-to surface mapping object
    const volSurfaceMapping vsm(aMesh);

    volScalarField Cvf
    (
        IOobject
        (
            "Cvf",
            runTime.timeName(),
            mesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar(dimless/dimLength, Zero)
    );

    vsm.mapToVolume(Cs, Cvf.boundaryFieldRef());
    Cvf.write();

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

    vsm.mapToVolume(Us, U.boundaryFieldRef());
    U.write();
