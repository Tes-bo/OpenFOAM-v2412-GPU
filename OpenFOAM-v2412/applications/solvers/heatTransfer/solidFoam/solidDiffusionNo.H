scalar DiNum = -GREAT;

{
    tmp<volScalarField> magKappa;
    if (thermo.isotropic())
    {
        magKappa = thermo.kappa();
    }
    else
    {
        magKappa = mag(thermo.Kappa());
    }

    tmp<volScalarField> tcp = thermo.Cp();
    const volScalarField& cp = tcp();

    tmp<volScalarField> trho = thermo.rho();
    const volScalarField& rho = trho();

    DiNum = max
    (
        solidRegionDiffNo
        (
            mesh,
            runTime,
            rho*cp,
            magKappa()
        ),
        DiNum
    );
}
