Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

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
    thermo.rho()
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

#include "compressibleCreatePhi.H"

pressureControl pressureControl(p, rho, pimple.dict(), false);

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

mesh.setFluxRequired(p.name());

#include "createDpdt.H"

#include "createK.H"

//- Overset specific

// Add solver-specific interpolations
{
    wordHashSet& nonInt =
        const_cast<wordHashSet&>(Stencil::New(mesh).nonInterpolatedFields());

    nonInt.insert("HbyA");
    nonInt.insert("grad(p)");
    nonInt.insert("surfaceIntegrate(phi)");
    nonInt.insert("surfaceIntegrate(phiHbyA)");
    nonInt.insert("cellMask");
    nonInt.insert("cellDisplacement");
    nonInt.insert("interpolatedCells");
    nonInt.insert("cellInterpolationWeight");
}

// Mask field for zeroing out contributions on hole cells
#include "createCellMask.H"
// Create bool field with interpolated cells
#include "createInterpolatedCells.H"

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);
