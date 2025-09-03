    Info<< "Reading boxTurbDict\n" << endl;

    IOdictionary boxTurbDict
    (
        IOobject
        (
            "boxTurbDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    scalar Ea(boxTurbDict.get<scalar>("Ea"));

    scalar k0(boxTurbDict.get<scalar>("k0"));
