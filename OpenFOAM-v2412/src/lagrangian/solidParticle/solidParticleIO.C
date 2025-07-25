/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2022 OpenCFD Ltd.
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

#include "solidParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::solidParticle::sizeofFields
(
    sizeof(solidParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::solidParticle::solidParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> d_ >> U_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, &d_);
            readRawScalar(is, U_.data(), vector::nComponents);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&d_), sizeofFields);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::solidParticle::readFields(Cloud<solidParticle>& c)
{
    const bool readOnProc = c.size();

    particle::readFields(c);

    IOField<scalar> d(c.newIOobject("d", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, d);

    IOField<vector> U(c.newIOobject("U", IOobject::MUST_READ), readOnProc);
    c.checkFieldIOobject(c, U);

    label i = 0;
    for (solidParticle& p : c)
    {
        p.d_ = d[i];
        p.U_ = U[i];
        ++i;
    }
}


void Foam::solidParticle::writeFields(const Cloud<solidParticle>& c)
{
    particle::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<scalar> d(c.newIOobject("d", IOobject::NO_READ), np);
    IOField<vector> U(c.newIOobject("U", IOobject::NO_READ), np);

    label i = 0;
    for (const solidParticle& p : c)
    {
        d[i] = p.d_;
        U[i] = p.U_;
        ++i;
    }

    d.write(writeOnProc);
    U.write(writeOnProc);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const solidParticle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.d_
            << token::SPACE << p.U_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.d_),
            solidParticle::sizeofFields
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
