{
    for (multiphaseInter::phaseModel& phase : fluid.phases())
    {
        PtrList<volScalarField>& Y = phase.Y();

        if (!Y.empty())
        {
            //- Su phase source terms
            PtrList<volScalarField::Internal> Sus(Y.size());
            //- Sp phase source terms
            PtrList<volScalarField::Internal> Sps(Y.size());

            forAll(Sus, i)
            {
                Sus.set
                (
                    i,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            "Su" + phase.name(),
                            mesh.time().timeName(),
                            mesh
                        ),
                        mesh,
                        dimensionedScalar(dimless/dimTime, Zero)
                    )
                );
                Sps.set
                (
                    i,
                    new volScalarField::Internal
                    (
                        IOobject
                        (
                            "Sp" + phase.name(),
                            mesh.time().timeName(),
                            mesh
                        ),
                        mesh,
                        dimensionedScalar(dimless/dimTime, Zero)
                    )
                );
            }

            forAll(Y, i)
            {
                // Calculate mass exchange for species consistent with
                // alpha's source terms.
                fluid.massSpeciesTransfer(phase, Sus[i], Sps[i], Y[i].name());
            }
            phase.solveYi(Sus, Sps);
        }
    }
}
