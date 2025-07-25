Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo(rhoThermo::New(mesh));
rhoThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    thermo.rho()
);

volScalarField& p = thermo.p();

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

autoPtr<compressible::turbulenceModel> turbulence;

IOobject turbulencePropertiesHeader
(
    "turbulenceProperties",
    runTime.constant(),
    mesh,
    IOobject::MUST_READ_IF_MODIFIED,
    IOobject::NO_WRITE,
    IOobject::NO_REGISTER
);

if (turbulencePropertiesHeader.typeHeaderOk<IOdictionary>(false))
{
    Info<< "Creating turbulence model\n" << endl;

    turbulence =
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        );
}

#include "createDpdt.H"

#include "createK.H"

#include "createRadiationModel.H"
#include "createFvOptions.H"
