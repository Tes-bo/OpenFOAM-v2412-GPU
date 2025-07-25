PtrList<uniformDimensionedScalarField> cumulativeContErrIO(fluidRegions.size());
forAll(cumulativeContErrIO, i)
{
    const fvMesh& mesh = fluidRegions[i];

    cumulativeContErrIO.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "cumulativeContErr",
                runTime.timeName(),
                "uniform",
                mesh.thisDb(),
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            dimensionedScalar(word::null, dimless, Zero)
        )
    );
}

UPtrList<scalar> cumulativeContErr(cumulativeContErrIO.size());
forAll(cumulativeContErrIO, i)
{
    cumulativeContErr.set(i, &cumulativeContErrIO[i].value());
}
