gmvFile << "tracers " << particles.size() << nl;
for (const passiveParticle& p : particles)
{
    gmvFile << p.position().x() << ' ';
}
gmvFile << nl;

for (const passiveParticle& p : particles)
{
    gmvFile << p.position().y() << ' ';
}
gmvFile << nl;

for (const passiveParticle& p : particles)
{
    gmvFile << p.position().z() << ' ';
}
gmvFile << nl;

for (const word& name : lagrangianScalarNames)
{
    IOField<scalar> fld
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (fld.size())
    {
        gmvFile << name  << nl;

        for (const scalar& val : fld)
        {
            gmvFile << val << token::SPACE;
        }
        gmvFile << nl;
    }


}

for (const word& name : lagrangianVectorNames)
{
    IOField<vector> fld
    (
        IOobject
        (
            name,
            runTime.timeName(),
            cloud::prefix,
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    if (fld.size())
    {
        gmvFile << name + "x" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].x() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "y" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].y() << token::SPACE;
        }
        gmvFile << nl;

        gmvFile << name + "z" << nl;

        forAll(fld, n)
        {
            gmvFile << fld[n].z() << token::SPACE;
        }
        gmvFile << nl;
    }

}


gmvFile << "endtrace"<< nl;
