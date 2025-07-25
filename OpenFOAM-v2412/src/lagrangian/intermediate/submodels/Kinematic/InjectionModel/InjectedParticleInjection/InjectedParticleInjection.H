/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2019 OpenCFD Ltd.
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
    Foam::InjectedParticleInjection

Description
    Replays an set of particle data based on an injectedParticleCloud, using
    the assumption of one particle per parcel.

Usage
    \verbatim
    model1
    {
        type            injectedParticleInjection;
        SOI             0;
        massTotal       0; // Place holder only
        parcelBasisType fixed;
        nParticle       1; // 1 particle per parcel
        cloud           eulerianParticleCloud;
        positionOffset  (-0.025 2 -0.025);
    }
    \endverbatim

SourceFiles
    InjectedParticleInjection.C

See also
    Foam::injectedParticle
    Foam::injectedParticleCloud
    Foam::functionObjects::extractEulerianParticles

\*---------------------------------------------------------------------------*/

#ifndef InjectedParticleInjection_H
#define InjectedParticleInjection_H

#include "InjectionModel.H"
#include "Switch.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class injectedParticleCloud;

/*---------------------------------------------------------------------------*\
                 Class InjectedParticleInjection Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class InjectedParticleInjection
:
    public InjectionModel<CloudType>
{
protected:

    // Protected data

        //- Name of cloud used to seed the new particles
        const word cloudName_;

        //- List of cell label per injector
        labelList injectorCells_;

        //- List of tetFace label per injector
        labelList injectorTetFaces_;

        //- List of tetPt label per injector
        labelList injectorTetPts_;

        //- List of injection time per particle [s]
        scalarList time_;

        //- List of position per particle [m]
        vectorList position_;

        //- Position offset to apply to input positions
        vector positionOffset_;

        //- List of diameter per particle [m]
        scalarList diameter_;

        //- List of velocity per particle [m/s]
        vectorList U_;

        //- List of volume per particle [m3]
        scalarList volume_;

        //- Flag to suppress errors if particle injection site is out-of-bounds
        Switch ignoreOutOfBounds_;

        //- Index of current particle
        label currentParticlei_;


        // Protected Member Functions

            //- Initialise injectors
            void initialise();


public:

    //- Runtime type information
    TypeName("injectedParticleInjection");


    // Constructors

        //- Construct from dictionary
        InjectedParticleInjection
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        InjectedParticleInjection
        (
            const InjectedParticleInjection<CloudType>& im
        );

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new InjectedParticleInjection<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~InjectedParticleInjection() = default;


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
                const label parceli,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            //- Set the parcel properties
            virtual void setProperties
            (
                const label parceli,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parceli);


        // I-O

            //- Write injection info
            void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "InjectedParticleInjection.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
