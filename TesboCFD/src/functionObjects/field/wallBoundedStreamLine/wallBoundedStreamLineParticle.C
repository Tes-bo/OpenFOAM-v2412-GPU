/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "wallBoundedStreamLineParticle.H"
#include "vectorFieldIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::vector Foam::wallBoundedStreamLineParticle::interpolateFields
(
    const trackingData& td,
    const point& position,
    const label celli,
    const label facei
)
{
    if (celli < 0)
    {
        FatalErrorInFunction
            << "Invalid cell (-1)" << abort(FatalError);
    }

    bool foundU = false;
    vector U(Zero);

    // If current position is different
    if
    (
        sampledPositions_.empty()
     || magSqr(sampledPositions_.last() - position) > Foam::sqr(SMALL)
    )
    {
        // Store new location
        sampledPositions_.append(position);

        // Scalar fields
        sampledScalars_.resize(td.vsInterp_.size());
        forAll(td.vsInterp_, i)
        {
            sampledScalars_[i].append
            (
                td.vsInterp_[i].interpolate(position, celli, facei)
            );
        }

        // Vector fields
        sampledVectors_.resize(td.vvInterp_.size());
        forAll(td.vvInterp_, i)
        {
            sampledVectors_[i].append
            (
                td.vvInterp_[i].interpolate(position, celli, facei)
            );

            if (td.vvInterp_.get(i) == &(td.UInterp_))
            {
                foundU = true;
                U = sampledVectors_[i].last();
            }
        }
    }

    if (!foundU)
    {
        U = td.UInterp_.interpolate(position, celli, facei);
    }

    return U;
}


Foam::vector Foam::wallBoundedStreamLineParticle::sample
(
    trackingData& td
)
{
    vector U = interpolateFields(td, localPosition_, cell(), face());

    const scalar magU = mag(U);

    if (magU < SMALL)
    {
        // Stagnant particle. Might as well stop
        lifeTime_ = 0;
        return vector::zero;
    }

    if (!trackForward_)
    {
        U = -U;
    }

    return U/magU;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    const point& position,
    const label celli,
    const label tetFacei,
    const label tetPti,
    const label meshEdgeStart,
    const label diagEdge,
    const bool trackForward,
    const label lifeTime
)
:
    wallBoundedParticle
    (
        mesh,
        position,
        celli,
        tetFacei,
        tetPti,
        meshEdgeStart,
        diagEdge
    ),
    trackForward_(trackForward),
    lifeTime_(lifeTime)
{}


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const polyMesh& mesh,
    Istream& is,
    bool readFields,
    bool newFormat
)
:
    wallBoundedParticle(mesh, is, readFields, newFormat)
{
    if (readFields)
    {
        List<scalarList> sampledScalars;
        List<vectorList> sampledVectors;

        is  >> trackForward_ >> lifeTime_
            >> sampledPositions_ >> sampledScalars >> sampledVectors;

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


Foam::wallBoundedStreamLineParticle::wallBoundedStreamLineParticle
(
    const wallBoundedStreamLineParticle& p
)
:
    wallBoundedParticle(p),
    trackForward_(p.trackForward_),
    lifeTime_(p.lifeTime_),
    sampledPositions_(p.sampledPositions_),
    sampledScalars_(p.sampledScalars_),
    sampledVectors_(p.sampledVectors_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::wallBoundedStreamLineParticle::readFields
(
    Cloud<wallBoundedStreamLineParticle>& c
)
{
    if (!c.size())
    {
        return;
    }

    wallBoundedParticle::readFields(c);

    IOField<label> lifeTime
    (
        c.newIOobject("lifeTime", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, lifeTime);

    vectorFieldIOField sampledPositions
    (
        c.newIOobject("sampledPositions", IOobject::MUST_READ)
    );
    c.checkFieldIOobject(c, sampledPositions);

    label i = 0;
    for (wallBoundedStreamLineParticle& p : c)
    {
        p.lifeTime_ = lifeTime[i];
        p.sampledPositions_.transfer(sampledPositions[i]);
        ++i;
    }
}


void Foam::wallBoundedStreamLineParticle::writeFields
(
    const Cloud<wallBoundedStreamLineParticle>& c
)
{
    wallBoundedParticle::writeFields(c);

    const label np = c.size();

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
    for (const wallBoundedStreamLineParticle& p : c)
    {
        lifeTime[i] = p.lifeTime_;
        sampledPositions[i] = p.sampledPositions_;
        ++i;
    }

    lifeTime.write();
    sampledPositions.write();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const wallBoundedStreamLineParticle& p
)
{
    os  << static_cast<const wallBoundedParticle&>(p)
        << token::SPACE << p.trackForward_
        << token::SPACE << p.lifeTime_
        << token::SPACE << p.sampledPositions_
        << token::SPACE << p.sampledScalars_
        << token::SPACE << p.sampledVectors_;

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
