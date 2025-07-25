{
    // Keep standard formulation on domain boundaries to ensure compatibility
    // with existing boundary conditions
    const Foam::FieldField<Foam::fvsPatchField, scalar> rhorAUf_orig
    (
        fvc::interpolate(rho.boundaryField()*rAU.boundaryField())
    );

    const Foam::FieldField<Foam::fvsPatchField, vector> rhoHbyA_orig
    (
        fvc::interpolate(rho.boundaryField()*HbyA.boundaryField())
    );

    surfaceScalarField::Boundary& rhorAUfbf = rhorAUf.boundaryFieldRef();
    surfaceVectorField::Boundary& rhoHbyAfbf = rhoHbyAf.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if (!U.boundaryField()[patchi].coupled())
        {
            rhorAUfbf[patchi] = rhorAUf_orig[patchi];
            rhoHbyAfbf[patchi] = rhoHbyA_orig[patchi];
        }
    }
}
