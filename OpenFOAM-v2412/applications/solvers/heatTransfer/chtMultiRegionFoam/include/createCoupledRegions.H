bool coupled = false;
{
    fvSolution solutionDict(runTime);
    solutionDict.readIfPresent("coupledEnergyField", coupled);
}

autoPtr<fvMatrix<scalar>> fvMatrixAssemblyPtr;

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(solidRegions, i)
{
    solidThermo& thermo = thermos[i];
    const auto& bpsi = thermo.T().boundaryField();

    forAll(bpsi, patchI)
    {
        if (bpsi[patchI].useImplicit())
        {
            coupled = true;
        }
    }
}

forAll(fluidRegions, i)
{
    const rhoThermo& thermo = refCast<const rhoThermo>(thermoFluid[i]);
    if (coupled)
    {
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}

if (coupled && !fvMatrixAssemblyPtr)
{
    forAll(solidRegions, i)
    {
        solidThermo& thermo = thermos[i];
        Info<< "Create fvMatrixAssembly." << endl;
        fvMatrixAssemblyPtr.reset
        (
            new fvMatrix<scalar>
            (
                thermo.he(),
                dimEnergy/dimTime
            )
        );
        break;
    }
}
