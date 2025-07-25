/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017-2022 OpenCFD Ltd.
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

#include "DTRMParticle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::DTRMParticle::propertyList_ =
    Foam::DTRMParticle::propertyList();

const std::size_t Foam::DTRMParticle::sizeofFields_
(
    sizeof(DTRMParticle) - sizeof(particle)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::DTRMParticle::DTRMParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    particle(mesh, is, readFields, newFormat),
    p0_(Zero),
    p1_(Zero),
    I0_(0),
    I_(0),
    dA_(0),
    transmissiveId_(-1)
{
    if (readFields)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is >> p0_ >> p1_ >> I0_ >> I_ >> dA_ >> transmissiveId_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, p0_.data(), vector::nComponents);
            readRawScalar(is, p1_.data(), vector::nComponents);
            readRawScalar(is, &I0_);
            readRawScalar(is, &I_);
            readRawScalar(is, &dA_);
            readRawLabel(is, &transmissiveId_);

            is.endRawRead();
        }
        else
        {
            is.read(reinterpret_cast<char*>(&p0_), sizeofFields_);
        }
    }

    is.check(FUNCTION_NAME);
}


void Foam::DTRMParticle::writeProperties
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

    writeProp("p0", p0_);
    writeProp("p1", p1_);
    writeProp("I0", I0_);
    writeProp("I", I_);
    writeProp("dA", dA_);
    writeProp("transmissiveId", transmissiveId_);

    #undef writeProp
}


Foam::Ostream& Foam::operator<<(Ostream& os, const DTRMParticle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.p0_
            << token::SPACE << p.p1_
            << token::SPACE << p.I0_
            << token::SPACE << p.I_
            << token::SPACE << p.dA_
            << token::SPACE << p.transmissiveId_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.p0_),
            DTRMParticle::sizeofFields_
        );
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);

    return os;
}


// ************************************************************************* //
