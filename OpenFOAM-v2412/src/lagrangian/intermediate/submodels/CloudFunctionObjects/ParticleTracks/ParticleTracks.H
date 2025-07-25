/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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
    Foam::ParticleTracks

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Records particle state (all variables) on each call to postFace

SourceFiles
    ParticleTracks.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleTracks_H
#define ParticleTracks_H

#include "CloudFunctionObject.H"
#include "labelPairHashes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleTracks Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleTracks
:
    public CloudFunctionObject<CloudType>
{
        //- Convenience typedef for parcel type
        typedef typename CloudType::parcelType parcelType;


    // Private Data

        //- Number of face-hit intervals between storing parcel data
        label trackInterval_;

        //- Maximum number of particles to store per track
        label maxSamples_;

        //- Should data should be reset/cleared on writing?
        bool resetOnWrite_;

        //- Table of number of times a particle has hit a face
        labelPairLookup faceHitCounter_;

        //- Pointer to the cloud storage
        autoPtr<Cloud<parcelType>> cloudPtr_;


protected:

    // Protected Member Functions

        //- Write post-processing info
        void write();


public:

    //- Runtime type information
    TypeName("particleTracks");


    // Constructors

        //- Construct from dictionary
        ParticleTracks
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        ParticleTracks(const ParticleTracks<CloudType>& ppm);

        //- Construct and return a clone
        virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
        {
            return autoPtr<CloudFunctionObject<CloudType>>
            (
                new ParticleTracks<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ParticleTracks() = default;


    // Member Functions

    // Access

        //- Return const access to the track interval
        inline label trackInterval() const;

        //- Return const access to the max samples
        inline label maxSamples() const;

        //- Should data be reset on write?
        inline bool resetOnWrite() const;

        //- Return the table of number of times a particle has hit a face
        inline const labelPairLookup& faceHitCounter() const;

        //- Return const access to the cloud
        inline const Cloud<parcelType>& cloud() const;


    // Evaluation

        //- Pre-evolve hook
        virtual void preEvolve
        (
            const typename parcelType::trackingData& td
        );

        //- Post-face hook
        virtual bool postFace
        (
            const parcelType& p,
            const typename parcelType::trackingData& td
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ParticleTracksI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ParticleTracks.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
