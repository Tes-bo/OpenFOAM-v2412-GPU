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
const scalar nu2 =
    dimensionedScalar
    (
        "nu",
        dimViscosity,
        transportProperties
    ).value();

Info<< "Reading viscoelastic properties\n" << endl;
IOdictionary turbulenceProperties
(
    IOobject
    (
        "turbulenceProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
const dictionary& MaxwellCoeffs =
    turbulenceProperties.subDict("laminar").subDict("MaxwellCoeffs");
const scalar nu1 = MaxwellCoeffs.get<scalar>("nuM");
const scalar lambda = MaxwellCoeffs.get<scalar>("lambda");

const scalar rho = 1;

Info<< "Reading pressure gradient\n" << endl;
IOdictionary fvOptions
(
    IOobject
    (
        "fvOptions",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

scalar K(0);

// Get x() component from U source
{
    const dictionary& momSource = fvOptions.subDict("momentumSource");

    if (momSource.findDict("sources"))
    {
        K = momSource.subDict("sources")
            .get<Tuple2<vector, scalar>>("U").first().x();
    }
    else
    {
        K = momSource.subDict("injectionRateSuSp")
            .get<Tuple2<vector, scalar>>("U").first().x();
    }
}

dictionary probes(IFstream(runTime.system()/"probes")());
const point location = pointField(probes.lookup("probeLocations"))[0];
