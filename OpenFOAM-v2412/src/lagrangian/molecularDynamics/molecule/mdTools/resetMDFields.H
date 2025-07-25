if (runTime.writeTime())
{
    allSpeciesN_RU = List<scalarField>
    (
        molecules.potential().nIds(),
        scalarField(mesh.nCells(), Zero)
    );

    allSpeciesM_RU = List<scalarField>
    (
        molecules.potential().nIds(),
        scalarField(mesh.nCells(), Zero)
    );

    allSpeciesVelocitySum_RU = List<vectorField>
    (
        molecules.potential().nIds(),
        vectorField(mesh.nCells(), Zero)
    );

    allSpeciesVelocityMagSquaredSum_RU = List<scalarField>
    (
        molecules.potential().nIds(),
        scalarField(mesh.nCells(), Zero)
    );
}
