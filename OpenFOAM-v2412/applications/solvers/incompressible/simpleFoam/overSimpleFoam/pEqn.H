{
    volScalarField rAU(1.0/UEqn.A());
    surfaceScalarField rAUf("rAUf", faceMask*fvc::interpolate(rAU));

    volVectorField HbyA("HbyA", U);
    HbyA = constrainHbyA(cellMask*rAU*UEqn.H(), U, p);

    tUEqn.clear();

    surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));

    adjustPhi(phiHbyA, U, p);

    if (adjustFringe)
    {
        oversetAdjustPhi(phiHbyA, U);
    }

    // Non-orthogonal pressure corrector loop
    while (simple.correctNonOrthogonal())
    {
        fvScalarMatrix pEqn
        (
            fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
        );

        pEqn.setReference(pRefCell, pRefValue);

        pEqn.solve();

        if (simple.finalNonOrthogonalIter())
        {
            phi = phiHbyA - pEqn.flux();
        }
    }

    #include "continuityErrs.H"

    // Explicitly relax pressure for momentum corrector
    p.relax();

    // Momentum corrector
    volVectorField gradP(fvc::grad(p));
    //mesh.interpolate(gradP);

    U = HbyA - rAU*cellMask*gradP;
    U.correctBoundaryConditions();
    fvOptions.correct(U);
}
