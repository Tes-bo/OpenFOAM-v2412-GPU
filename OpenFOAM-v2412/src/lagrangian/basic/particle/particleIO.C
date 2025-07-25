/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "particle.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::string Foam::particle::propertyList_  = Foam::particle::propertyList();

const std::size_t Foam::particle::sizeofPosition
(
    offsetof(particle, facei_) - offsetof(particle, coordinates_)
);

const std::size_t Foam::particle::sizeofFields
(
    sizeof(particle) - offsetof(particle, coordinates_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::particle::particle
(
    const polyMesh& mesh,
    Istream& is,
    const bool readFields,
    const bool newFormat,
    const bool doLocate
)
:
    mesh_(mesh),
    coordinates_(),
    celli_(-1),
    tetFacei_(-1),
    tetPti_(-1),
    facei_(-1),
    stepFraction_(0.0),
    origProc_(Pstream::myProcNo()),
    origId_(-1)
{
    point position;
    readData(is, position, readFields, newFormat, doLocate);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::particle::readData
(
    Istream& is,
    point& position,
    const bool readFields,
    const bool newFormat,
    const bool doLocate
)
{
    if (newFormat)
    {
        if (is.format() == IOstreamOption::ASCII)
        {
            is  >> coordinates_ >> celli_ >> tetFacei_ >> tetPti_;
            if (readFields)
            {
                is  >> facei_ >> stepFraction_ >> origProc_ >> origId_;
            }
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, coordinates_.data(), barycentric::nComponents);
            readRawLabel(is, &celli_);
            readRawLabel(is, &tetFacei_);
            readRawLabel(is, &tetPti_);

            if (readFields)
            {
                readRawLabel(is, &facei_);
                readRawScalar(is, &stepFraction_);
                readRawLabel(is, &origProc_);
                readRawLabel(is, &origId_);
            }

            is.endRawRead();
        }
        else
        {
            if (readFields)
            {
                is.read(reinterpret_cast<char*>(&coordinates_), sizeofFields);
            }
            else
            {
                is.read(reinterpret_cast<char*>(&coordinates_), sizeofPosition);
            }
        }
    }
    else
    {
        positionsCompat1706 p;

        if (is.format() == IOstreamOption::ASCII)
        {
            is >> p.position >> p.celli;

            if (readFields)
            {
                is  >> p.facei
                    >> p.stepFraction
                    >> p.tetFacei
                    >> p.tetPti
                    >> p.origProc
                    >> p.origId;
            }
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size

            is.beginRawRead();

            readRawScalar(is, p.position.data(), vector::nComponents);
            readRawLabel(is, &p.celli);

            if (readFields)
            {
                readRawLabel(is, &p.facei);
                readRawScalar(is, &p.stepFraction);
                readRawLabel(is, &p.tetFacei);
                readRawLabel(is, &p.tetPti);
                readRawLabel(is, &p.origProc);
                readRawLabel(is, &p.origId);
            }

            is.endRawRead();
        }
        else
        {
            if (readFields)
            {
                // Read whole struct
                const size_t s =
                (
                    sizeof(positionsCompat1706)
                  - offsetof(positionsCompat1706, position)
                );
                is.read(reinterpret_cast<char*>(&p.position), s);
            }
            else
            {
                // Read only position and cell
                const size_t s =
                (
                    offsetof(positionsCompat1706, facei)
                  - offsetof(positionsCompat1706, position)
                );
                is.read(reinterpret_cast<char*>(&p.position), s);
            }
        }

        if (readFields)
        {
            // Note: other position-based properties are set using locate(...)
            stepFraction_ = p.stepFraction;
            origProc_ = p.origProc;
            origId_ = p.origId;
        }

        // Preserve read position
        position = p.position;

        if (doLocate)
        {
            locate
            (
                p.position,
                nullptr,
                p.celli,
                false,
                "Particle initialised with a location outside of the mesh."
            );
        }
    }

    // Check state of Istream
    is.check(FUNCTION_NAME);
}


void Foam::particle::writeProperties
(
    Ostream& os,
    const wordRes& filters,
    const word& delim,
    const bool namesOnly
) const
{
    #undef  writeProp
    #define writeProp(Name, Value)                                            \
        particle::writeProperty(os, Name, Value, namesOnly, delim, filters)

    writeProp("coordinates", coordinates_);
    writeProp("position", position());
    writeProp("celli", celli_);
    writeProp("tetFacei", tetFacei_);
    writeProp("tetPti", tetPti_);
    writeProp("facei", facei_);
    writeProp("stepFraction", stepFraction_);
    writeProp("origProc", origProc_);
    writeProp("origId", origId_);

    #undef writeProp
}


void Foam::particle::writeCoordinates(Ostream& os) const
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << coordinates_
            << token::SPACE << celli_
            << token::SPACE << tetFacei_
            << token::SPACE << tetPti_;
    }
    else
    {
        os.write(reinterpret_cast<const char*>(&coordinates_), sizeofPosition);
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


void Foam::particle::writePosition(Ostream& os) const
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << position() << token::SPACE << celli_;
    }
    else
    {
        positionsCompat1706 p;

        const size_t s =
        (
            offsetof(positionsCompat1706, facei)
          - offsetof(positionsCompat1706, position)
        );

        p.position = position();
        p.celli = celli_;

        os.write(reinterpret_cast<const char*>(&p.position), s);
    }

    // Check state of Ostream
    os.check(FUNCTION_NAME);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const particle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << p.coordinates_
            << token::SPACE << p.celli_
            << token::SPACE << p.tetFacei_
            << token::SPACE << p.tetPti_
            << token::SPACE << p.facei_
            << token::SPACE << p.stepFraction_
            << token::SPACE << p.origProc_
            << token::SPACE << p.origId_;
    }
    else
    {
        os.write
        (
            reinterpret_cast<const char*>(&p.coordinates_),
            particle::sizeofFields
        );
    }

    return os;
}


// ************************************************************************* //
