if (pimple.nCorrPIMPLE() > 1)
{
    if (!pimple.firstIter())
    {
        // Resetting alpha1 to value before advection in first PIMPLE
        // iteration.
        alpha1 = alpha1.oldTime();
    }
}
tmp<surfaceScalarField> talphaPhi1(alphaPhi10);

if (nAlphaSubCycles > 1)
{
    dimensionedScalar totalDeltaT = runTime.deltaT();

    talphaPhi1.reset
    (
        new surfaceScalarField
        (
            mesh.newIOobject("alphaPhi1"),
            mesh,
            dimensionedScalar(alphaPhi10.dimensions(), Zero)
        )
    );

    surfaceScalarField rhoPhiSum
    (
        mesh.newIOobject("rhoPhiSum"),
        mesh,
        dimensionedScalar(rhoPhi.dimensions(), Zero)
    );

    for
    (
        subCycle<volScalarField> alphaSubCycle(alpha1, nAlphaSubCycles);
        !(++alphaSubCycle).end();
    )
    {
        #include "alphaEqn.H"
        talphaPhi1.ref() += (runTime.deltaT()/totalDeltaT)*alphaPhi10;
        rhoPhiSum += (runTime.deltaT()/totalDeltaT)*rhoPhi;
    }

    rhoPhi = rhoPhiSum;
}
else
{
    #include "alphaEqn.H"
}

rho == alpha1*rho1 + alpha2*rho2;

const surfaceScalarField& alphaPhi1 = talphaPhi1();
surfaceScalarField alphaPhi2("alphaPhi2", phi - alphaPhi1);

volScalarField::Internal contErr
(
    (
        fvc::ddt(rho) + fvc::div(rhoPhi)
      - (fvOptions(alpha1, mixture.thermo1().rho())&rho1)
      - (fvOptions(alpha2, mixture.thermo2().rho())&rho2)
    )()
);
