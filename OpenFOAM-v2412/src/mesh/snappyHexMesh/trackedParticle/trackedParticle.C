/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "trackedParticle.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const std::size_t Foam::trackedParticle::sizeofFields_
(
    sizeof(trackedParticle) - offsetof(trackedParticle, start_)
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::trackedParticle::trackedParticle
(
    const polyMesh& mesh,
    const barycentric& coordinates,
    const label celli,
    const label tetFacei,
    const label tetPtI,
    const point& end,
    const label level,
    const label i,
    const label j,
    const label k
)
:
    particle(mesh, coordinates, celli, tetFacei, tetPtI),
    start_(position()),
    end_(end),
    level_(level),
    i_(i),
    j_(j),
    k_(k)
{}


Foam::trackedParticle::trackedParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const point& end,
    const label level,
    const label i,
    const label j,
    const label k
)
:
    particle(mesh, position, celli),
    start_(this->position()),
    end_(end),
    level_(level),
    i_(i),
    j_(j),
    k_(k)
{}


Foam::trackedParticle::trackedParticle
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
            is >> start_ >> end_ >> level_ >> i_ >> j_ >> k_;
        }
        else if (!is.checkLabelSize<>() || !is.checkScalarSize<>())
        {
            // Non-native label or scalar size
            is.beginRawRead();

            readRawScalar(is, start_.data(), vector::nComponents);
            readRawScalar(is, end_.data(), vector::nComponents);
            readRawLabel(is, &level_);
            readRawLabel(is, &i_);
            readRawLabel(is, &j_);
            readRawLabel(is, &k_);

            is.endRawRead();
        }
        else
        {
            is.read
            (
                reinterpret_cast<char*>(&start_),
                sizeofFields_
            );
        }
    }

    is.check(FUNCTION_NAME);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::trackedParticle::move
(
    Cloud<trackedParticle>& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    td.switchProcessor = false;

    scalar tEnd = (1.0 - stepFraction())*trackTime;

    if (tEnd <= SMALL && onBoundaryFace())
    {
        // This is a hack to handle particles reaching their endpoint
        // on a processor boundary. If the endpoint is on a processor face
        // it currently gets transferred backwards and forwards infinitely.

        // Remove the particle
        td.keepParticle = false;
    }
    else
    {
        td.keepParticle = true;

        while (td.keepParticle && !td.switchProcessor && stepFraction() < 1)
        {
            // mark visited cell with max level.
            td.maxLevel_[cell()] = max(td.maxLevel_[cell()], level_);

            const scalar f = 1 - stepFraction();
            const vector s = end_ - start_;
            trackToAndHitFace(f*s, f, cloud, td);
        }
    }

    return td.keepParticle;
}


bool Foam::trackedParticle::hitPatch(Cloud<trackedParticle>&, trackingData&)
{
    return false;
}


void Foam::trackedParticle::hitWedgePatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPlanePatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitSymmetryPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicAMIPatch
(
    Cloud<trackedParticle>&,
    trackingData& td,
    const vector& direction
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitCyclicACMIPatch
(
    Cloud<trackedParticle>&,
    trackingData& td,
    const vector&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::hitProcessorPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Move to different processor
    td.switchProcessor = true;
}


void Foam::trackedParticle::hitWallPatch
(
    Cloud<trackedParticle>&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::trackedParticle::correctAfterParallelTransfer
(
    const label patchi,
    trackingData& td
)
{
    particle::correctAfterParallelTransfer(patchi, td);

    label edgeI = k();
    if (edgeI != -1)
    {
        label featI = i();

        // Mark edge we're currently on (was set on sending processor but not
        // receiving sender)
        td.featureEdgeVisited_[featI].set(edgeI);
    }
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const trackedParticle& p)
{
    if (os.format() == IOstreamOption::ASCII)
    {
        os  << static_cast<const particle&>(p)
            << token::SPACE << p.start_
            << token::SPACE << p.end_
            << token::SPACE << p.level_
            << token::SPACE << p.i_
            << token::SPACE << p.j_
            << token::SPACE << p.k_;
    }
    else
    {
        os  << static_cast<const particle&>(p);
        os.write
        (
            reinterpret_cast<const char*>(&p.start_),
            trackedParticle::sizeofFields_
        );
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
