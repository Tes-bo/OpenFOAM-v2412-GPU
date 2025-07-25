    Info<< "Reading field h" << endl;
    areaScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );


    Info<< "Reading field Us" << endl;
    areaVectorField Us
    (
        IOobject
        (
            "Us",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
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
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(Us) & aMesh.Le()
    );


    edgeScalarField phi2s
    (
        IOobject
        (
            "phi2s",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        fac::interpolate(h*Us) & aMesh.Le()
    );


    const areaVectorField& Ns = aMesh.faceAreaNormals();
    areaVectorField Gs(g - Ns*(Ns & g));
    areaScalarField Gn(mag(g - Gs));

    // Mass source
    areaScalarField Sm
    (
        IOobject
        (
            "Sm",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimLength/dimTime, Zero)
    );

    // Mass sink
    areaScalarField Sd
    (
        IOobject
        (
            "Sd",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar(dimLength/dimTime, Zero)
    );

    areaVectorField Ug
    (
        IOobject
        (
            "Sg",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedVector(dimVelocity, Zero)
    );


    // Surface pressure
    areaScalarField ps
    (
        IOobject
        (
            "ps",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rhol*Gn*h - sigma*fac::laplacian(h)
    );

    // Friction factor
    areaScalarField dotProduct
    (
        aMesh.faceAreaNormals() & (g/mag(g))
    );

    Info<< "View factor: min = " << min(dotProduct.internalField())
        << " max = " << max(dotProduct.internalField()) << endl;

    areaScalarField manningField
    (
        IOobject
        (
            "manningField",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    areaScalarField frictionFactor
    (
        IOobject
        (
            "frictionFactor",
            runTime.timeName(),
            aMesh.thisDb(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh,
        dimensionedScalar("one", dimless, 0.01)
    );

    aMesh.setFluxRequired("h");
