    Info<< "Reading field T\n" << endl;

    volScalarField T
    (
        IOobject
        (
            "T",
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


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT("DT", dimViscosity, transportProperties);

    bool oversetPatchErrOutput =
        simple.dict().getOrDefault("oversetPatchErrOutput", false);

    // Dummy phi for oversetPatchErrOutput
    tmp<surfaceScalarField> tdummyPhi;
    if (oversetPatchErrOutput)
    {
        tdummyPhi = tmp<surfaceScalarField>::New
        (
            IOobject
            (
                "dummyPhi",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, Zero)
        );
    }
