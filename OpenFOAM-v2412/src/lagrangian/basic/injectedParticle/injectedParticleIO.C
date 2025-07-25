/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2022 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "injectedParticle.H"
#include "IOstreams.H"
#include "IOField.H"
#include "Cloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::injectedParticle::propertyList_ =
    Foam::injectedParticle::propertyList();


const std::size_t Foam::injectedParticle::sizeofFields
(
    // Does not include position_
    sizeof(injectedParticle) - offsetof(injectedParticle, tag_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::injectedParticle::injectedParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, false), // force to read old positions file
    position_(Zero),
    tag_(-1),
    soi_(0.0),
    d_(0.0),
    U_(Zero)
{
    if (readFields)
    {
        // After the base particle class has read the fields from file and
        // constructed the necessary barycentric coordinates we can update the
        // particle position on this mesh
        position_ = particle::position();

        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> tag_ >> soi_ >> d_ >> U_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size
            is.beginRawRead();

            readRawLabel(is, &tag_);
            readRawScalar(is, &soi_);
            readRawScalar(is, &d_);
            readRawScalar(is, U_.data(), vector::nComponents);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&tag_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::injectedParticle::readFields(Cloud<injectedParticle>& c)
{
    if (!c.size())
    {
        return;
    }

    // Note: not reading local position_ - defer to base particle class

    particle::readFields(c);

    IOField<label> tag(c.newIOobject("tag", IOobject::MUST_READ));
    c.checkFieldIOobject(c, tag);

    IOField<scalar> soi(c.newIOobject("soi", IOobject::MUST_READ));
    c.checkFieldIOobject(c, soi);

    IOField<scalar> d(c.newIOobject("d", IOobject::MUST_READ));
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.newIOobject("U", IOobject::MUST_READ));
    c.checkFieldIOobject(c, U);

    label i = 0;
    for (injectedParticle& p : c)
    {
        p.tag_ = tag[i];
        p.soi_ = soi[i];
        p.d_ = d[i];
        p.U_ = U[i];

        ++i;
    }
}


void Foam::injectedParticle::writeFields(const Cloud<injectedParticle>& c)
{
    // Force writing positions instead of coordinates
    const bool oldWriteCoordinates = particle::writeLagrangianCoordinates;
    const bool oldWritePositions = particle::writeLagrangianPositions;

    particle::writeLagrangianCoordinates = false;
    particle::writeLagrangianPositions = true;

    particle::writeFields(c);

    // Note: not writing local position_ - defer to base particle class

    const label np = c.size();

    IOField<label> tag(c.newIOobject("tag", IOobject::NO_READ), np);
    IOField<scalar> soi(c.newIOobject("soi", IOobject::NO_READ), np);
    IOField<scalar> d(c.newIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.newIOobject("U", IOobject::NO_READ), np);

    label i = 0;

    for (const injectedParticle& p : c)
    {
        tag[i] = p.tag();
        soi[i] = p.soi();
        d[i] = p.d();
        U[i] = p.U();

        ++i;
    }

    tag.write();
    soi.write();
    d.write();
    U.write();

    // Restore
    particle::writeLagrangianCoordinates = oldWriteCoordinates;
    particle::writeLagrangianPositions = oldWritePositions;
}


void Foam::injectedParticle::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    particle::writeProperties(os, filters, delim, namesOnly);

    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        particle::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("tag", tag_);
    writeProp("soi", soi_);
    writeProp("d", d_);
    writeProp("U", U_);

    #undef  writeProp
}


void Foam::injectedParticle::readObjects
(
    Cloud<injectedParticle>& c,
    const objectRegistry& obr
)
{
    particle::readObjects(c, obr);

    if (!c.size()) return;

    const auto& tag = cloud::lookupIOField<label>("tag", obr);
    const auto& soi = cloud::lookupIOField<scalar>("soi", obr);
    const auto& d = cloud::lookupIOField<scalar>("d", obr);
    const auto& U = cloud::lookupIOField<vector>("U", obr);

    label i = 0;

    for (injectedParticle& p : c)
    {
         p.tag() = tag[i];
         p.soi() = soi[i];
         p.d() = d[i];
         p.U() = U[i];

        ++i;
    }
}


void Foam::injectedParticle::writeObjects
(
    const Cloud<injectedParticle>& c,
    objectRegistry& obr
)
{
    // Always writes "position", not "coordinates"
    particle::writeObjects(c, obr);

    const label np = c.size();

    auto& tag = cloud::createIOField<label>("tag", np, obr);
    auto& soi = cloud::createIOField<scalar>("soi", np, obr);
    auto& d = cloud::createIOField<scalar>("d", np, obr);
    auto& U = cloud::createIOField<vector>("U", np, obr);

    label i = 0;

    for (const injectedParticle& p : c)
    {
        tag[i] = p.tag();
        soi[i] = p.soi();
        d[i] = p.d();
        U[i] = p.U();

        ++i;
    }
}


void Foam::injectedParticle::writePosition(Ostream& os) const
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << position_ << token::SPACE << cell();
    }
    else
    {
        positionsCompat1706 p;

        const size_t s =
        (
            offsetof(positionsCompat1706, facei)
          - offsetof(positionsCompat1706, position));

        p.position = position_;
        p.celli = cell();

        os.write(reinterpret_cast<const char*>(&p.position), s);
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const injectedParticle& p
)
{
    // Note: not writing local position_ - defer to base particle class

    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.tag()
            << token::SPACE << p.soi()
            << token::SPACE << p.d()
            << token::SPACE << p.U();
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.tag_),
            injectedParticle::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
