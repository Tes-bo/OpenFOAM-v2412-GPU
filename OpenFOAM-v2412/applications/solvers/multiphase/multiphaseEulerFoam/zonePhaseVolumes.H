{
    const scalarField& V = mesh.V();

    forAll(mesh.cellZones(), czi)
    {
        const labelList& cellLabels = mesh.cellZones()[czi];

        for (const volScalarField& alpha : fluid.phases())
        {
            scalar phaseVolume = 0;

            for (const label celli : cellLabels)
            {
                phaseVolume += alpha[celli]*V[celli];
            }

            reduce(phaseVolume, sumOp<scalar>());

            Info<< alpha.name()
                << " phase volume in zone " << mesh.cellZones()[czi].name()
                << " = " << phaseVolume*1e6 << " ml " << endl;
        }
    }
}
