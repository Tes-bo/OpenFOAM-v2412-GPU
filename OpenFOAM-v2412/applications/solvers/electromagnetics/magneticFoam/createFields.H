    Info<< "Reading field psi\n" << endl;
    volScalarField psi
    (
        IOobject
        (
            "psi",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading transportProperties\n" << endl;

    IOdictionary transportProperties
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    List<magnet> magnets(transportProperties.lookup("magnets"));

    surfaceScalarField murf
    (
        IOobject
        (
            "murf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("one", dimless, 1.0)
    );

    surfaceScalarField Mrf
    (
        IOobject
        (
            "Mrf",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar(dimensionSet(0, 1, 0, 0, 0, 1, 0), Zero)
    );
    Mrf.setOriented(true);

    forAll(magnets, i)
    {
        label magnetZonei = mesh.faceZones().findZoneID(magnets[i].name());

        if (magnetZonei == -1)
        {
            FatalIOErrorInFunction(transportProperties)
                << "Cannot find faceZone for magnet " << magnets[i].name()
                << exit(FatalIOError);
        }

        const labelList& faces = mesh.faceZones()[magnetZonei];

        const scalar muri = magnets[i].mur();
        const scalar Mri = magnets[i].Mr().value();
        const vector& orientationi = magnets[i].orientation();

        const surfaceVectorField& Sf = mesh.Sf();

        for (const label facei : faces)
        {
            murf[facei] = muri;
            Mrf[facei] = Mri*(orientationi & Sf[facei]);
        }
    }
