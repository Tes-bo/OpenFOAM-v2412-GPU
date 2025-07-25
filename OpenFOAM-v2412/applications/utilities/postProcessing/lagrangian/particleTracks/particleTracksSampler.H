/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

Class
    Foam::particleTracksSampler

Description
    Helper class when generating particle tracks.
    The interface is fairly rudimentary.

SourceFiles
    particleTracksSamplerTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_particleTracksSampler_H
#define Foam_particleTracksSampler_H

#include "globalIndex.H"
#include "fieldTypes.H"
#include "IOField.H"
#include "labelField.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class objectRegistry;

/*---------------------------------------------------------------------------*\
                    Class particleTracksSampler Declaration
\*---------------------------------------------------------------------------*/

class particleTracksSampler
{
    // Private Data

        //- The sampling interval
        label stride_ = 1;

        //- The number of tracks
        label nTracks_ = 1;

        //- Upper limit for number of track positions
        label maxPositions_ = 10000;

        //- The bin for sampled parcels.
        labelField trackIds_;

        //- The addressing for gathering cloud fields etc
        globalIndex cloudGather_;

        //- The originating parcel addressing
        globalIndex origParcelAddr_;


public:

    // Data Members

        //- The originating parcel ids
        labelField origParcelIds_;

        //- The originating processor ids
        labelField origProcIds_;


    // Member Functions

        //- The original parcel addressing
        const globalIndex& parcelAddr() const noexcept
        {
            return origParcelAddr_;
        }

        //- Total number of particles
        label nParticle() const
        {
            return origParcelAddr_.totalSize();
        }

        //- Number of tracks to generate
        label nTracks() const noexcept
        {
            return nTracks_;
        }

        //- Define the orig parcel mappings
        void reset(const labelUList& origParcelCounts)
        {
            origParcelAddr_.reset(origParcelCounts);
            origParcelIds_.clear();
            origProcIds_.clear();
            trackIds_.clear();
        }

        //- Set the sampling stride, upper limits
        //  \return Number of tracks to generate
        label setSampleRate
        (
            const label sampleFreq,
            const label maxPositions,
            const label maxTracks = -1
        )
        {
            // numTracks = numParticles/stride

            stride_ = max(1, sampleFreq);
            maxPositions_ = maxPositions;
            nTracks_ = (origParcelAddr_.totalSize()/stride_);

            if (maxTracks > 0)
            {
                nTracks_ = min(nTracks_, maxTracks);
            }

            return nTracks_;
        }

        void resetCloud(const label localCloudSize)
        {
            cloudGather_.reset(globalIndex::gatherOnly{}, localCloudSize);
            origParcelIds_.resize_nocopy(localCloudSize);
            origProcIds_.resize_nocopy(localCloudSize);
        }

        template<class Type>
        void gatherInplace(List<Type>& fld) const
        {
            cloudGather_.gatherInplace(fld);
        }

        template<class Type>
        void createTrackField
        (
            const UList<Type>& values,
            List<DynamicList<Type>>& trackValues
        ) const;


        template<class Type>
        label setTrackFields
        (
            const objectRegistry& obr,
            HashTable<List<DynamicList<Type>>>& fieldTable
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "particleTracksSamplerTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
