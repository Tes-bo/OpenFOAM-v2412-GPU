Info<< "Constructing face momentum equations" << endl;

PtrList<fvVectorMatrix> UEqns(phases.size());

{
    fluid.momentumTransfer(); // !!! Update coefficients shouldn't be necessary
                              //     This should be done on demand

    autoPtr<phaseSystem::momentumTransferTable>
        momentumTransferPtr(fluid.momentumTransferf());

    phaseSystem::momentumTransferTable&
        momentumTransfer(momentumTransferPtr());

    forAll(fluid.movingPhases(), movingPhasei)
    {
        phaseModel& phase = fluid.movingPhases()[movingPhasei];

        const volScalarField& alpha = phase;
        const tmp<volScalarField> trho = phase.rho();
        volVectorField& U = phase.URef();

        const auto& rho = trho();

        UEqns.set
        (
            phase.index(),
            new fvVectorMatrix
            (
                phase.UfEqn()
             ==
               *momentumTransfer[phase.name()]
              + fvOptions(alpha, rho, U)
            )
        );

        UEqns[phase.index()].relax();
        fvOptions.constrain(UEqns[phase.index()]);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}
