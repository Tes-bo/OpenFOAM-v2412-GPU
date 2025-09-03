/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2023 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TrackCloudType>
bool Foam::wallBoundedStreamLineParticle::move
(
    TrackCloudType& cloud,
    trackingData& td,
    const scalar trackTime
)
{
    auto& p = static_cast<typename TrackCloudType::particleType&>(*this);

    // Check position is inside tet
    //checkInside();

    td.switchProcessor = false;
    td.keepParticle = true;

    scalar tEnd = (1.0 - stepFraction())*trackTime;
    scalar maxDt = mesh().bounds().mag();

    while (td.keepParticle && !td.switchProcessor && lifeTime_ > 0)
    {
        // Set the lagrangian time-step
        scalar dt = maxDt;

        --lifeTime_;

        // Get sampled velocity and fields. Store if position changed.
        vector U = p.sample(td);

        // !user parameter!
        if (dt < SMALL)
        {
            // Force removal
            lifeTime_ = 0;
            break;
        }


        if (td.trackLength_ < GREAT)
        {
            dt = td.trackLength_;
        }


        scalar fraction = trackToEdge(cloud, td, localPosition_ + dt*U);
        dt *= fraction;

        tEnd -= dt;
        stepFraction() = 1.0 - tEnd/trackTime;


        if (tEnd <= ROOTVSMALL)
        {
            // Force removal
            lifeTime_ = 0;
        }
    }


    if (!td.keepParticle || lifeTime_ == 0)
    {
        if (lifeTime_ == 0)
        {
            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle :"
                    << " Removing stagnant particle:"
                    << localPosition_
                    << " sampled positions:" << sampledPositions_.size()
                    << endl;
            }
            td.keepParticle = false;
        }
        else
        {
            // Normal exit. Store last position and fields
            sample(td);

            if (debug)
            {
                Pout<< "wallBoundedStreamLineParticle : Removing particle:"
                    << localPosition_
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


// ************************************************************************* //
