for (phaseModel& phase : fluid.phases())
{
    phase.DDtU() =
        fvc::ddt(phase.U())
      + fvc::div(phase.phi(), phase.U())
      - fvc::div(phase.phi())*phase.U();

    MRF.addAcceleration(phase.U(), phase.DDtU());
}
