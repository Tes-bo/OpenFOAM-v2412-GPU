    Info<< "Reading field p_rgh\n" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    Info<< "Reading field U\n" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Note: construct T to be around before the thermos. The thermos will
    //       not update T.
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


    Info<< "Calculating field g.h\n" << endl;
    #include "readGravitationalAcceleration.H"
    #include "readhRef.H"
    #include "gh.H"

    volScalarField  p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );


    Info<< "Creating multiphaseSystem\n" << endl;
    autoPtr<multiphaseInter::multiphaseSystem> fluidPtr =
        multiphaseInter::multiphaseSystem::New(mesh);

    multiphaseInter::multiphaseSystem& fluid = fluidPtr();

    if (!fluid.incompressible())
    {
        FatalError << "One or more phases are not incompressible. " << nl
            << "This is a incompressible solver." << abort(FatalError);
    }


    // Need to store rho for ddt(rho, U)
    volScalarField rho
    (
        IOobject
        (
            "rho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        fluid.rho()
    );
    rho.oldTime();

    // Update p using fluid.rho()
    p = p_rgh + rho*gh;

    label pRefCell = 0;
    scalar pRefValue = 0.0;
    setRefCell
    (
        p,
        p_rgh,
        pimple.dict(),
        pRefCell,
        pRefValue
    );

    if (p_rgh.needReference())
    {
        p += dimensionedScalar
        (
            "p",
            p.dimensions(),
            pRefValue - getRefCellValue(p, pRefCell)
        );
        p_rgh = p - rho*gh;
    }

    // Mass flux
    surfaceScalarField& rhoPhi = fluid.rhoPhi();

    // Construct incompressible turbulence model
    autoPtr
    <
        CompressibleTurbulenceModel<multiphaseInter::multiphaseSystem>
    > turbulence
    (
        CompressibleTurbulenceModel<multiphaseInter::multiphaseSystem>::New
        (
            rho,
            U,
            rhoPhi,
            fluid
        )
    );

    fluid.setTurbulence(turbulence());

    // Creating radiation model
    autoPtr<radiation::radiationModel> radiation
    (
        radiation::radiationModel::New(T)
    );

    Info<< "Calculating field rhoCp\n" << endl;
    volScalarField rhoCp
    (
        IOobject
        (
            "rhoCp",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fluid.rho()*fluid.Cp()
    );
    rhoCp.oldTime();


    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(U));
