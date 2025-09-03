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

gmvFile << "U" << nl;
for (const passiveParticle& p : particles)
{
    gmvFile << p.velocity().x() << " ";
}
gmvFile << nl;

gmvFile << "V" << nl;
for (const passiveParticle& p : particles)
{
    gmvFile << p.velocity().y() << " ";
}
gmvFile << nl;

gmvFile << "W" << nl;
for (const passiveParticle& p : particles)
{
{
    gmvFile << p.velocity().z() << " ";
}
gmvFile << nl;

gmvFile << "Diam" << nl;
for (const passiveParticle& p : particles)
{
    gmvFile << p.d() << " ";
}

gmvFile << "endtrace"<< nl;
