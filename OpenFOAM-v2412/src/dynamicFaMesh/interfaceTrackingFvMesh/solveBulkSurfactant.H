{
    volScalarField& C = bulkSurfactantConcentration();

    const dimensionedScalar& D = surfactant().bulkDiffusion();

    scalar ka = surfactant().adsorptionCoeff().value();
    scalar kb = surfactant().desorptionCoeff().value();
    scalar CsInf = surfactant().saturatedConc().value();

    const scalarField& Cs =
        surfactantConcentration().internalField();

    // const scalarField& Cfs = C.boundaryField()[fsPatchIndex()];

    if
    (
        C.boundaryField()[fsPatchIndex()].type()
     == fixedGradientFvPatchScalarField::typeName
    )
    {
        fixedGradientFvPatchScalarField& fsC =
            refCast<fixedGradientFvPatchScalarField>
            (
                C.boundaryFieldRef()[fsPatchIndex()]
            );

        fsC.gradient() = (ka*kb*Cs - ka*(CsInf-Cs)*fsC)/D.value();
    }
    else
    {
        FatalErrorInFunction
            << "Bulk concentration boundary condition "
            << "at the free-surface patch is not "
            << fixedGradientFvPatchScalarField::typeName
            << exit(FatalError);
    }

    fvScalarMatrix CEqn
    (
        fvm::ddt(C)
      + fvm::div(phi(), C, "div(phi,C)")
      - fvm::laplacian(D, C, "laplacian(D,C)")
    );

    CEqn.solve();
}
