{
    alphav = clamp((rho - rholSat)/(rhovSat - rholSat), zero_one{});
    alphal = 1.0 - alphav;

    Info<< "max-min alphav: " << max(alphav).value()
        << " " << min(alphav).value() << endl;

    psiModel->correct();

    //Info<< "min a: " << 1.0/sqrt(max(psi)).value() << endl;
}
