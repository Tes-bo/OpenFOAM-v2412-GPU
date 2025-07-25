#include "MRFCorrectBCs.H"

PtrList<fvVectorMatrix> UEqns(fluid.phases().size());
autoPtr<multiphaseSystem::dragCoeffFields> dragCoeffs(fluid.dragCoeffs());

label phasei = 0;
for (phaseModel& phase : fluid.phases())
{
    const volScalarField& alpha = phase;
    volVectorField& U = phase.U();

    volScalarField nuEff(turbulence->nut() + phase.nu());

    UEqns.set
    (
        phasei,
        new fvVectorMatrix
        (
            fvm::ddt(alpha, U)
          + fvm::div(phase.alphaPhi(), U)

          + (alpha/phase.rho())*fluid.Cvm(phase)*
            (
                fvm::ddt(U)
              + fvm::div(phase.phi(), U)
              - fvm::Sp(fvc::div(phase.phi()), U)
            )

          - fvm::laplacian(alpha*nuEff, U)
          - fvc::div
            (
                alpha*(nuEff*dev(T(fvc::grad(U))) /*- ((2.0/3.0)*I)*k*/),
                "div(Rc)"
            )
         ==
        //- fvm::Sp(fluid.dragCoeff(phase, dragCoeffs())/phase.rho(), U)
        //- (alpha*phase.rho())*fluid.lift(phase)
          //+
            (alpha/phase.rho())*fluid.Svm(phase)
          - fvm::Sp
            (
                slamDampCoeff
               *max
                (
                    mag(U()) - maxSlamVelocity,
                    dimensionedScalar("U0", dimVelocity, Zero)
                )
               /cbrt(mesh.V()),
                U
            )
        )
    );
    MRF.addAcceleration
    (
        alpha*(1 + (1/phase.rho())*fluid.Cvm(phase)),
        UEqns[phasei]
    );
    //UEqns[phasei].relax();

    ++phasei;
}
