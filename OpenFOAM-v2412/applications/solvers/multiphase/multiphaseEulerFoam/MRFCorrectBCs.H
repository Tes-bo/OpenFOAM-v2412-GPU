    for (phaseModel& phase : fluid.phases())
    {
        MRF.correctBoundaryVelocity(phase.U());
    }

    MRF.correctBoundaryVelocity(U);
