IOobject dpdtHeader
(
    "dpdt",
    runTime.timeName(),
    mesh,
    IOobject::NO_READ,
    IOobject::NO_WRITE
);

if (mesh.dynamic())
{
    Info<< "Creating field dpdt for moving meshes\n" << endl;

    // Note
    // - set to READ_IF_PRESENT and AUTO_WRITE to simplify dpdt correction
    //   by meshPhi

    dpdtHeader.readOpt(IOobject::READ_IF_PRESENT);
    dpdtHeader.writeOpt(IOobject::AUTO_WRITE);
}
else
{
    Info<< "Creating field dpdt\n" << endl;
}

volScalarField dpdt(dpdtHeader, fvc::ddt(p));

if (!thermo.dpdt())
{
    dpdt == dimensionedScalar(dpdt.dimensions(), Zero);
    dpdt.writeOpt(IOobject::NO_WRITE);
}
