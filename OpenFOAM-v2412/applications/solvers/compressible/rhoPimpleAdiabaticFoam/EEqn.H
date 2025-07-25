{
    volScalarField& he = thermo.he();

    const tmp<volScalarField>& tCp = thermo.Cp();
    const tmp<volScalarField>& tCv = thermo.Cv();

    const volScalarField& Cp = tCp();
    const volScalarField& Cv = tCv();
    const scalar gamma = max(Cp/Cv).value();

    if (mag(gamma - min(Cp/Cv).value()) > VSMALL)
    {
        notImplemented("gamma not constant in space");
    }

    const dictionary& thermoDict = thermo.subDict("mixture");

    const dictionary& eosDict = thermoDict.subDict("equationOfState");

    bool local = eosDict.getOrDefault("local", false);

    // Evolve T as:
    //
    // T_1 = T_0 \frac{p}{p_0}^{\frac{\gamma - 1}{\gamma}}

    if (!local)
    {
        const scalar T0 = eosDict.get<scalar>("T0");
        const scalar p0 = eosDict.get<scalar>("p0");

        he = thermo.he(p, pow(p/p0, (gamma - scalar(1))/gamma)*T0);
    }
    else
    {
        const volScalarField& T0 = T.oldTime();
        const volScalarField& p0 = p.oldTime();

        he = thermo.he(p, pow(p/p0, (gamma - scalar(1))/gamma)*T0);
    }

    thermo.correct();

    psi = 1.0/((Cp - Cv)*T);

    rho = thermo.rho();
    rho.relax();

    rho.writeMinMax(Info);
}
