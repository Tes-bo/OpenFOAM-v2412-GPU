// Initialise fluid field pointer lists
PtrList<twoPhaseSystem> phaseSystemFluid(fluidRegions.size());

PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());

PtrList<volScalarField> p_rghFluid(fluidRegions.size());

PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());
List<bool> frozenFlowFluid(fluidRegions.size(), false);

List<label> pRefCellFluid(fluidRegions.size());
List<scalar> pRefValueFluid(fluidRegions.size());
PtrList<dimensionedScalar> pMinFluid(fluidRegions.size());

const uniformDimensionedVectorField& g = meshObjects::gravity::New(runTime);

PtrList<pimpleControl> pimpleFluid(fluidRegions.size());


//Debug Fields
/*
PtrList<volScalarField> faceRegimesFluid(fluidRegions.size());
PtrList<volScalarField> qcFluid(fluidRegions.size());
PtrList<volScalarField> qFilmFluid(fluidRegions.size());
PtrList<volScalarField> htcFilmBoilingFluid(fluidRegions.size());
PtrList<volScalarField> qtbFluid(fluidRegions.size());
PtrList<volScalarField> qSubCoolFluid(fluidRegions.size());
PtrList<volScalarField> CHFtotalFluid(fluidRegions.size());
PtrList<volScalarField> TdnbFluid(fluidRegions.size());
PtrList<volScalarField> phiFluid(fluidRegions.size());

forAll(fluidRegions, i)
{
    faceRegimesFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "faceRegimes",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)

        )
    );
    qcFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qc",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    qFilmFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qFilm",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    htcFilmBoilingFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "htcFilmBoiling",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    qtbFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qtb",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    qSubCoolFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "qSubCool",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
    CHFtotalFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "CHFtotal",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
     TdnbFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Tdnb",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
     phiFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "phiTb",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar(dimless, Zero)
        )
    );
}
*/

forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    pimpleFluid.set
    (
        i,
        new pimpleControl(fluidRegions[i])
    );

    Info<< "    Adding to phaseSystemFluid\n" << endl;
    phaseSystemFluid.set(i, twoPhaseSystem::New(fluidRegions[i]).ptr());

    Info<< "    Adding hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(word::null, dimLength, Zero)
        )
    );

    Info<< "    Adding ghRef\n" << endl;
    dimensionedScalar ghRef
    (
        mag(g.value()) > SMALL
      ? g & (cmptMag(g.value())/mag(g.value()))*hRefFluid[i]
      : dimensionedScalar("ghRef", g.dimensions()*dimLength, 0)
    );

    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (g & fluidRegions[i].C()) - ghRef
        )
    );

    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (g & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding p_rghFluid\n" << endl;
    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Correcting p_rghFluid\n" << endl;

    // Force p_rgh to be consistent with p
    p_rghFluid[i] =
        phaseSystemFluid[i].phase1().thermo().p()
     -  phaseSystemFluid[i].phase1().thermo().rho()*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Correcting initialMassFluid\n" << endl;
    initialMassFluid[i] =
        fvc::domainIntegrate(phaseSystemFluid[i].rho()).value();

    const dictionary& pimpleDict =
        fluidRegions[i].solutionDict().subDict("PIMPLE");

    pimpleDict.readIfPresent("frozenFlow", frozenFlowFluid[i]);

    pRefCellFluid[i] = -1;
    pRefValueFluid[i] = 0.0;

    Info<< "    Setting reference\n" << endl;
    if (p_rghFluid[i].needReference())
    {
        setRefCell
        (
            phaseSystemFluid[i].phase1().thermoRef().p(),
            p_rghFluid[i],
            pimpleDict,
            pRefCellFluid[i],
            pRefValueFluid[i]
        );
    }

    pMinFluid.set
    (
        i,
        new dimensionedScalar
        (
            "pMin",
            dimPressure,
            phaseSystemFluid[i]
        )
    );
}
