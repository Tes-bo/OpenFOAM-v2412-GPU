{
    volScalarField rAU(1.0/UEqn.A());
    volVectorField HbyA("HbyA", U);
    HbyA = rAU*UEqn.H();


    // Define coefficients and pseudo-velocities for RCM interpolation
    // M[U] = AU - H = -grad(p)
    // U = H/A - 1/A grad(p)
    // H/A = U + 1/A grad(p)
    surfaceScalarField rhorAUf
    (
        "rhorAUf",
        fvc::interpolate(rho)/fvc::interpolate(UEqn.A())
    );

    surfaceVectorField rhoHbyAf
    (
        "rhoHbyAf",
        fvc::interpolate(rho)*fvc::interpolate(U)
      + rhorAUf*fvc::interpolate(fvc::grad(p))
    );

    #include "resetBoundaries.H"

    if (pimple.nCorrPISO() <= 1)
    {
        tUEqn.clear();
    }

    if (pimple.transonic())
    {
         FatalError
             << "\nTransonic option not available for " << args.executable()
             << exit(FatalError);
    }
    else
    {
        // Rhie & Chow interpolation (part 1)
        surfaceScalarField phiHbyA
        (
            "phiHbyA",
            (
                (rhoHbyAf & mesh.Sf())
              + rhorAUf*fvc::interpolate(rho)*fvc::ddtCorr(U, phiByRho)
              + fvc::interpolate(rho)
              * fvc::alphaCorr(U, phiByRho, pimple.finalInnerIter())
            )
        );

        MRF.makeRelative(fvc::interpolate(rho), phiHbyA);

        // Non-orthogonal pressure corrector loop
        while (pimple.correctNonOrthogonal())
        {
            // Pressure corrector
            fvScalarMatrix pEqn
            (
                fvm::ddt(psi, p)
              + fvc::div(phiHbyA)
              - fvm::laplacian(rhorAUf, p)
              ==
                fvOptions(psi, p, rho.name())
            );

            pEqn.solve(p.select(pimple.finalInnerIter()));

            // Rhie & Chow interpolation (part 2)
            if (pimple.finalNonOrthogonalIter())
            {
                phi = phiHbyA + pEqn.flux();
            }
        }
    }

    phiByRho = phi/fvc::interpolate(rho);

    #include "rhoEqn.H"
    #include "compressibleContinuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    U = HbyA - rAU*fvc::grad(p);
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}

rho = thermo.rho();
