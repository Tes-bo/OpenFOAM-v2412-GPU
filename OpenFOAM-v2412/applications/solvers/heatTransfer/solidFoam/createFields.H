Info<< "Reading thermophysical properties\n" << endl;

autoPtr<solidThermo> pThermo(solidThermo::New(mesh));
solidThermo& thermo = pThermo();

tmp<volScalarField> trho = thermo.rho();

autoPtr<coordinateSystem> coordinatesPtr;
autoPtr<volSymmTensorField> taniAlpha;

if (!thermo.isotropic())
{
    Info<< "Adding coordinateSystem\n" << endl;
    coordinatesPtr = coordinateSystem::New
    (
        mesh,
        thermo,
        coordinateSystem::typeName
    );

    tmp<volVectorField> tkappaByCp = thermo.Kappa()/thermo.Cp();

    taniAlpha.reset
    (
        new volSymmTensorField
        (
            IOobject
            (
                "Anialpha",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedSymmTensor(tkappaByCp().dimensions(), Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    volSymmTensorField& aniAlpha = *taniAlpha;

    aniAlpha.primitiveFieldRef() =
        coordinatesPtr->transformPrincipal
        (
            mesh.cellCentres(),
            tkappaByCp()
        );
    aniAlpha.correctBoundaryConditions();
}


IOobject betavSolidIO
(
    "betavSolid",
    runTime.timeName(),
    mesh,
    IOobject::MUST_READ,
    IOobject::AUTO_WRITE
);

autoPtr<volScalarField> betavPtr;
if (betavSolidIO.typeHeaderOk<volScalarField>(true))
{
    betavPtr.reset
    (
        new volScalarField
        (
            betavSolidIO,
            mesh
        )
    );
}
else
{
    betavPtr.reset
    (
        new volScalarField
        (
            IOobject
            (
                "betavSolid",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("1", dimless, scalar(1))
        )
    );
}

// Consider mesh flux to correct for mesh deformation
bool meshFluxCorr(false);
if (mesh.solutionDict().found("SIMPLE"))
{
    meshFluxCorr =
        mesh.solutionDict().subDict("SIMPLE").getOrDefault<bool>
        (
            "meshFluxCorrection", false
        );
}
else if (mesh.solutionDict().found("PIMPLE"))
{
    meshFluxCorr =
        mesh.solutionDict().subDict("PIMPLE").getOrDefault<bool>
        (
            "meshFluxCorrection", false
        );
}

#include "createRadiationModel.H"
#include "createFvOptions.H"
