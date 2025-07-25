/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::InjectedParticleDistributionInjection

Description
    Interrogates an injectedParticleCloud to convert the raw particle data
    into a set of 'binned' injectors.

    The bins are set according to the particle \c tag property, from which:
    - diameters are converted into \c general distributions with a
      user-specified bin width
    - raw velocity and diameter data are resampled and stored to provide
      variations per injector

    The mass to inject can be set according to the raw input data mass total
    by using the \c applyDistributionMassTotal switch


Usage
    \verbatim
    model1
    {
        type            injectedParticleDistributionInjection;
        SOI             0;
        parcelBasisType mass;
        cloud           eulerianParticleCloud;
        positionOffset  (-0.025 2 -0.025);
        binWidth        0.1e-3;
        parcelsPerInjector 500;
        resampleSize    100;                    // optional
        applyDistributionMassTotal yes;

        // Placeholder only when using applyDistributionMassTotal
        massTotal       0;
    }
    \endverbatim

Note
    The each injector location is assumed to be operating under steady
    conditions, i.e. using a constant flow rate profile

SourceFiles
    InjectedParticleDistributionInjection.C

See also
    Foam::injectedParticle
    Foam::injectedParticleCloud
    Foam::functionObjects::extractEulerianParticles
    Foam::distributionModels::general

\*---------------------------------------------------------------------------*/

#ifndef InjectedParticleDistributionInjection_H
#define InjectedParticleDistributionInjection_H

#include "InjectionModel.H"
#include "general.H"
#include "Switch.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class InjectedParticleDistributionInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectedParticleDistributionInjection
:
    public InjectionModel<CloudType>
{
protected:

    // Protected data

        //- Name of cloud used to seed the new particles
        const word cloudName_;

        //- List of start time per injector
        scalarList startTime_;

        //- List of end time per injector
        scalarList endTime_;

        //- List of position per injector
        List<vectorList> position_;

        //- Position offset to apply to input positions
        vector positionOffset_;

        //- List of volume flow rate per injector [m3/s]
        scalarList volumeFlowRate_;

        //- List of parcel velocity per injector
        List<vectorList> U_;

        //- Bin width when generating particle distributions
        scalar binWidth_;

        //- List of size distribution model per injector
        PtrList<distributionModels::general> sizeDistribution_;

        //- Target number of parcels to inject per injector
        scalar parcelsPerInjector_;

        //- Resample size
        label resampleSize_;

        //- Flag to apply mass calculated from distribution instead of
        //  InjectionModel massTotal_
        bool applyDistributionMassTotal_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;

        //- Running total of number of parcels injected
        label nParcelsInjected_;

        //- Number of parcels injected in last step (local proc only)
        label nParcelsInjected0_;

        //- Current injector
        label currentInjectori_;

        //- Current sample
        label currentSamplei_;


    // Protected Member Functions

        //- Initialise injectors
        void initialise();


public:

    //- Runtime type information
    TypeName("injectedParticleDistributionInjection");


    // Constructors

        //- Construct from dictionary
        InjectedParticleDistributionInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        InjectedParticleDistributionInjection
        (
            const InjectedParticleDistributionInjection<CloudType>& im
        );

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new InjectedParticleDistributionInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InjectedParticleDistributionInjection();


    // Member Functions

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFaceI,
                label& tetPtI
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);


        // I-O

            //- Write injection info
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectedParticleDistributionInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
