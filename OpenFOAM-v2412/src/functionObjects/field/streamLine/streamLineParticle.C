/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "streamLineParticle.H"
#include "streamLineParticleCloud.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::streamLineParticle::interpolateFields
(
    const trackingData& td,
    const barycentric& tetCoords,
    const tetIndices& tetIs
)
{
    if (tetIs.cell() < 0)
    {
        FatalErrorInFunction
            << "Invalid cell (-1)" << abort(FatalError);
    }

    const point position
    (
        tetIs.barycentricToPoint(this->mesh(), tetCoords)
    );

    bool foundU = false;
    vector U(Zero);

    // If current position is different
    if
    (
        sampledPositions_.empty()
     || sampledPositions_.back().distSqr(position) > Foam::sqr(SMALL)
    )
    {
        // Store new location
        sampledPositions_.push_back(position);

        // Scalar fields
        sampledScalars_.resize(td.vsInterp_.size());
        forAll(td.vsInterp_, i)
        {
            sampledScalars_[i].push_back
            (
                td.vsInterp_[i].interpolate(tetCoords, tetIs, tetIs.face())
            );
        }

        // Vector fields
        sampledVectors_.resize(td.vvInterp_.size());
        forAll(td.vvInterp_, i)
        {
            sampledVectors_[i].push_back
            (
                td.vvInterp_[i].interpolate(tetCoords, tetIs, tetIs.face())
            );

            if (td.vvInterp_.get(i) == &(td.UInterp_))
            {
                foundU = true;
                U = sampledVectors_[i].back();
            }
        }
    }

    if (!foundU)
    {
        U = td.UInterp_.interpolate(tetCoords, tetIs, tetIs.face());
    }

    return U;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::streamLineParticle::streamLineParticle
(
    const polyMesh& mesh,
    const vector& position,
    const label celli,
    const bool trackForward,
    const label lifeTime
)
:
    particle(mesh, position, celli),
    trackForward_(trackForward),
    lifeTime_(lifeTime)
{}


Foam::streamLineParticle::streamLineParticle
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
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> trackForward_ >> lifeTime_
            >> sampledPositions_ >> sampledScalars
            >> sampledVectors;

        sampledScalars_.resize(sampledScalars.size());
        forAll(sampledScalars, i)
        {
            sampledScalars_[i].transfer(sampledScalars[i]);
        }
        sampledVectors_.resize(sampledVectors.size());
        forAll(sampledVectors, i)
        {
            sampledVectors_[i].transfer(sampledVectors[i]);
        }
    }

    is.check(FUNCTION_NAME);
}


Foam::streamLineParticle::streamLineParticle
(
    const streamLineParticle& p
)
:
    particle(p),
    trackForward_(p.trackForward_),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_),
    sampledVectors_(p.sampledVectors_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::streamLineParticle::move
(
    streamLineParticleCloud& cloud,
    trackingData& td,
    const scalar
)
{
    td.switchProcessor = false;
    td.keepParticle = true;

    const scalar maxDt = mesh().bounds().mag();

    while (td.keepParticle && !td.switchProcessor && lifeTime_ > 0)
    {
        // Set the lagrangian time-step
        scalar dt = maxDt;

        // Cross cell in steps:
        // - at subiter 0 calculate dt to cross cell in nSubCycle steps
        // - at the last subiter do all of the remaining track
        for (label subIter = 0; subIter < max(1, td.nSubCycle_); subIter++)
        {
            --lifeTime_;

            // Store current position, return sampled velocity
            vector U =
                interpolateFields(td, coordinates(), currentTetIndices());

            const scalar magU = mag(U);

            if (magU < SMALL)
            {
                // Stagnant particle. Might as well stop
                lifeTime_ = 0;
                break;
            }

            if (!trackForward_)
            {
                U = -U;
            }

            U /= magU;

            if (td.trackLength_ < GREAT)
            {
                // No sub-cycling. Track a set length on each step.
                dt = td.trackLength_;
            }
            else if (subIter == 0)
            {
                // Sub-cycling. Cross the cell in nSubCycle steps.
                particle copy(*this);
                copy.trackToFace(maxDt*U, 1);
                dt *= (copy.stepFraction() - stepFraction())/td.nSubCycle_;
            }
            else if (subIter == td.nSubCycle_ - 1)
            {
                // Sub-cycling. Track the whole cell on the last step.
                dt = maxDt;
            }

            trackToAndHitFace(dt*U, 0, cloud, td);

            if
            (
                onFace()
            || !td.keepParticle
            ||  td.switchProcessor
            ||  lifeTime_ == 0
            )
            {
                break;
            }
        }
    }

    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            // Failure exit. Particle stagnated or it's life ran out.
            if (debug)
            {
                Pout<< "streamLineParticle: Removing stagnant particle:"
                    << position() << " sampled positions:"
                    << sampledPositions_.size() << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            (void)interpolateFields(td, coordinates(), currentTetIndices());

            if (debug)
            {
                Pout<< "streamLineParticle: Removing particle:" << position()
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
        }

        // Transfer particle data into trackingData.
        {
            td.allPositions_.emplace_back().transfer(sampledPositions_);
        }

        forAll(sampledScalars_, i)
        {
            td.allScalars_[i].emplace_back().transfer(sampledScalars_[i]);
        }
        forAll(sampledVectors_, i)
        {
            td.allVectors_[i].emplace_back().transfer(sampledVectors_[i]);
        }
    }

    return td.keepParticle;
}


bool Foam::streamLineParticle::hitPatch(streamLineParticleCloud&, trackingData&)
{
    // Disable generic patch interaction
    return false;
}


void Foam::streamLineParticle::hitWedgePatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPlanePatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitSymmetryPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicAMIPatch
(
    streamLineParticleCloud&,
    trackingData& td,
    const vector&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitCyclicACMIPatch
(
    streamLineParticleCloud&,
    trackingData& td,
    const vector&
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::hitProcessorPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Switch particle
    td.switchProcessor = true;
}


void Foam::streamLineParticle::hitWallPatch
(
    streamLineParticleCloud&,
    trackingData& td
)
{
    // Remove particle
    td.keepParticle = false;
}


void Foam::streamLineParticle::readFields(Cloud<streamLineParticle>& c)
{
    const bool readOnProc = c.size();

    particle::readFields(c);

    IOField<label> lifeTime
    (
        c.newIOobject("lifeTime", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.newIOobject("sampledPositions", IOobject::MUST_READ),
        readOnProc
    );
    c.checkFieldIOobject(c, sampledPositions);

    label i = 0;
    for (streamLineParticle& p : c)
    {
        p.lifeTime_ = lifeTime[i];
        p.sampledPositions_.transfer(sampledPositions[i]);
        ++i;
    }
}


void Foam::streamLineParticle::writeFields(const Cloud<streamLineParticle>& c)
{
    particle::writeFields(c);

    const label np = c.size();
    const bool writeOnProc = c.size();

    IOField<label> lifeTime
    (
        c.newIOobject("lifeTime", IOobject::NO_READ),
        np
    );
    vectorFieldIOField sampledPositions
    (
        c.newIOobject("sampledPositions", IOobject::NO_READ),
        np
    );

    label i = 0;
    for (const streamLineParticle& p : c)
    {
        lifeTime[i] = p.lifeTime_;
        sampledPositions[i] = p.sampledPositions_;
        ++i;
    }

    lifeTime.write(writeOnProc);
    sampledPositions.write(writeOnProc);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const streamLineParticle& p)
{
    os  << static_cast<const particle&>(p)
        << token::SPACE << p.trackForward_
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
