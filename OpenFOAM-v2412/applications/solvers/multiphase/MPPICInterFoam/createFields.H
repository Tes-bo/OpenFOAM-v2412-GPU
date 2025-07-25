#include "createRDeltaT.H"

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

#include "createPhi.H"

Info<< "Reading transportProperties\n" << endl;
immiscibleIncompressibleTwoPhaseMixture mixture(U, phi);

volScalarField& alpha1(mixture.alpha1());
volScalarField& alpha2(mixture.alpha2());

const dimensionedScalar& rho1 = mixture.rho1();
const dimensionedScalar& rho2 = mixture.rho2();

// Need to store rho for ddt(rho, U)
volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alpha1*rho1 + alpha2*rho2
);
rho.oldTime();

// Need to store mu as incompressibleTwoPhaseMixture does not store it
volScalarField mu
(
    IOobject
    (
        "mu",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT
    ),
    mixture.mu()
);


// Mass flux
surfaceScalarField rhoPhi
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::interpolate(rho)*phi
);

#include "readGravitationalAcceleration.H"
#include "readhRef.H"
#include "gh.H"

volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    p_rgh + rho*gh
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell
(
    p,
    p_rgh,
    mesh.solutionDict().subDict("PIMPLE"),
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

mesh.setFluxRequired(p_rgh.name());
mesh.setFluxRequired(alpha1.name());

// alphac must be constructed before the cloud
// so that the drag-models can find it
volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar(dimless, Zero),
    fvPatchFieldBase::zeroGradientType()
);
alphac.oldTime();

volScalarField alphacRho(alphac*rho);
alphacRho.oldTime();

Info<< "Constructing kinematicCloud " << endl;
basicKinematicCloud kinematicCloud
(
    "kinematicCloud",
    rho,
    U,
    mu,
    g
);

// Particle fraction upper limit
scalar alphacMin
(
    1.0
  - (
        kinematicCloud.particleProperties().subDict("constantProperties")
       .get<scalar>("alphaMax")
    )
);

// Update alphac from the particle locations
alphac = max(1.0 - kinematicCloud.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphacf", fvc::interpolate(alphac));

// Phase mass flux
surfaceScalarField alphaRhoPhic("alphaRhoPhic", alphacf*rhoPhi);

// Volumetric phase flux
surfaceScalarField alphaPhic("alphaPhic", alphacf*phi);

autoPtr
<
    PhaseCompressibleTurbulenceModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    >
>turbulence
(
    PhaseCompressibleTurbulenceModel
    <
        immiscibleIncompressibleTwoPhaseMixture
    >::New
    (
        alphac,
        rho,
        U,
        alphaRhoPhic,
        rhoPhi,
        mixture
    )
);

#include "createMRF.H"
