if (adjustTimeStep)
{
    runTime.setDeltaT(min(dtChem, maxDeltaT));
    Info<< "deltaT = " <<  runTime.deltaTValue() << endl;
}
