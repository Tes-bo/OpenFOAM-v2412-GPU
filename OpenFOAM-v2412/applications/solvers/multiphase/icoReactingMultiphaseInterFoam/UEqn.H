fvVectorMatrix UEqn
(
    fvm::ddt(rho, U)
  + fvm::div(rhoPhi, U)
//- fvm::Sp(fvc::ddt(rho) + fvc::div(rhoPhi), U)
  + turbulence->divDevRhoReff(U)
  ==
    fvOptions(rho, U)
);

UEqn.relax();

fluid.addInterfacePorosity(UEqn);

if (pimple.momentumPredictor())
{
    solve
    (
        UEqn
      ==
        fvc::reconstruct
        (
            (
                fluid.surfaceTensionForce()
                - ghf*fvc::snGrad(rho)
                - fvc::snGrad(p_rgh)
            ) * mesh.magSf()
        )
    );

    fvOptions.correct(U);
    K = 0.5*magSqr(U);
}
