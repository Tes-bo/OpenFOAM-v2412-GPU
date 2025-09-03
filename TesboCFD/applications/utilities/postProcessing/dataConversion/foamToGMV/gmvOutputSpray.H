gmvFile << "tracers " << particles.size() << nl;
for (const passiveParticle& p : particles)
{
    gmvFile << p.position().x() << " ";
}
gmvFile << nl;

for (const passiveParticle& p : particles)
{
    gmvFile << p.position().y() << " ";
}
gmvFile << nl;

for (const passiveParticle& p : particles)
{
    gmvFile << p.position().z() << " ";
}
gmvFile << nl;

forAll(lagrangianScalarNames, i)
{
    word name = lagrangianScalarNames[i];

    IOField<scalar> s
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

    if (s.size())
    {
        gmvFile << name  << nl;

        forAll(s, n)
        {
            gmvFile << s[n] << token::SPACE;
        }
        gmvFile << nl;
    }

}

gmvFile << "endtrace"<< nl;
