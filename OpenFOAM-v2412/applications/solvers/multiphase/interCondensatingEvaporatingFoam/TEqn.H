{
    tmp<volScalarField> tcp(thermo->Cp());
    const volScalarField& cp = tcp();

    const dimensionedScalar Cp1 = thermo->Cp1();
    const dimensionedScalar Cp2 = thermo->Cp2();

    rhoCp = rho*cp;

    kappaEff = thermo->kappa() + rho*cp*turbulence->nut()/Prt;

    const surfaceScalarField rhoCpPhi
    (
        "rhoCpPhi",
        rhoPhi*(Cp1 - Cp2) + phi*rho2*Cp2
    );

    fvScalarMatrix TEqn
    (
         fvm::ddt(rhoCp, T)
       + fvm::div(rhoCpPhi, T)
       - fvm::Sp(fvc::ddt(rhoCp) + fvc::div(rhoCpPhi), T)
       - fvm::laplacian(kappaEff, T)
       + mixture->TSource()
    );


    TEqn.relax();
    TEqn.solve();

    Info<< "min/max(T) = " << min(T).value() << ", "
        << max(T).value() <<endl;
}
