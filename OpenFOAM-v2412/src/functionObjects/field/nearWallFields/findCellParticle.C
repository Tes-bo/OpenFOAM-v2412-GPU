/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2018-2022 OpenCFD Ltd.
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

#include "findCellParticle.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::findCellParticle::findCellParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const point& end,
    const label data
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPti),
    start_(position()),
    end_(end),
    data_(data)
{}


Foam::findCellParticle::findCellParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const point& end,
    const label data
)
:
    particle(mesh, position, celli),
    start_(this->position()),
    end_(end),
    data_(data)
{}


Foam::findCellParticle::findCellParticle
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
            is >> start_ >> end_ >> data_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size
            is.beginRawRead();

            readRawScalar(is, start_.data(), vector::nComponents);
            readRawScalar(is, end_.data(), vector::nComponents);
            readRawLabel(is, &data_);

            is.endRawRead();
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&start_),
                sizeof(start_) + sizeof(end_) + sizeof(data_)
            );
        }
    }

    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::findCellParticle::move
(
    Cloud<findCellParticle>& cloud,
    trackingData& td,
    const scalar maxTrackLen
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
    {
        const scalar f = 1 - stepFraction();
        trackToAndHitFace(f*(end_ - start_), f, cloud, td);
    }

    // Note: stepFraction is might not be exactly 1 so check for 1 or
    // slightly larger

    if (stepFraction() >= 1 || !td.keepParticle)
    {
        // Hit endpoint or patch. If patch hit could do fancy stuff but just
        // to use the patch point is good enough for now.
        td.cellToData()[cell()].append(data());
        td.cellToEnd()[cell()].append(position());
    }

    return td.keepParticle;
}


bool Foam::findCellParticle::hitPatch(Cloud<findCellParticle>&, trackingData&)
{
    return false;
}


void Foam::findCellParticle::hitWedgePatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPlanePatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitSymmetryPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitCyclicPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitCyclicAMIPatch
(
    Cloud<findCellParticle>&,
    trackingData& td,
    const vector&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitCyclicACMIPatch
(
    Cloud<findCellParticle>&,
    trackingData& td,
    const vector&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::findCellParticle::hitProcessorPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.switchProcessor = true;
}


void Foam::findCellParticle::hitWallPatch
(
    Cloud<findCellParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const findCellParticle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.start_
            << token::SPACE << p.end_
            << token::SPACE << p.data_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.start_),
            sizeof(p.start_) + sizeof(p.end_) + sizeof(p.data_)
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
