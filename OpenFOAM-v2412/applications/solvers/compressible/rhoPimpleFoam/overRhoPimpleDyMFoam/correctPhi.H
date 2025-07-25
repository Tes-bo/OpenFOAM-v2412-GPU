if (mesh.changing())
{
    volVectorField::Boundary& bfld = U.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].initEvaluate();
        }
    }

    surfaceScalarField::Boundary& phiBfld = phi.boundaryFieldRef();
    forAll(bfld, patchi)
    {
        if (bfld[patchi].fixesValue())
        {
            bfld[patchi].evaluate();

            phiBfld[patchi] =
                rho.boundaryField()[patchi]
              * (
                    bfld[patchi]
                  & mesh.Sf().boundaryField()[patchi]
                );
        }
    }
}
    // Initialize BCs list for pcorr to zero-gradient
    wordList pcorrTypes
    (
        p.boundaryField().size(),
        fvPatchFieldBase::zeroGradientType()
    );

    // Set BCs of pcorr to fixed-value for patches at which p is fixed
    forAll(p.boundaryField(), patchi)
    {
        if (p.boundaryField()[patchi].fixesValue())
        {
            pcorrTypes[patchi] = fixedValueFvPatchScalarField::typeName;
        }
    }

    volScalarField pcorr
    (
        IOobject
        (
            "pcorr",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(p.dimensions(), Zero),
        pcorrTypes
    );

    mesh.setFluxRequired(pcorr.name());

{
    dimensionedScalar rAUf("rAUf", dimTime, 1.0);

    while (pimple.correctNonOrthogonal())
    {
        fvScalarMatrix pcorrEqn
        (
            fvm::ddt(psi, pcorr)
          + fvc::div(phi)
          - fvm::laplacian(rAUf, pcorr)
         ==
            divrhoU()
        );

        pcorrEqn.solve(pcorr.select(pimple.finalInnerIter()));
        //Bypass virtual layer
        //mesh.fvMesh::solve(pcorrEqn, d);

        if (pimple.finalNonOrthogonalIter())
        {
            phi += pcorrEqn.flux();
        }
    }
}
