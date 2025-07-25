/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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
    Foam::KinematicSurfaceFilm

Group
    grpLagrangianIntermediateSurfaceFilmSubModels

Description
    Kinematic parcel surface film model.

    Responsible for:
    - injecting parcels from the film model into the cloud, e.g. for dripping
    - parcel interaction with the film, e.g absorb, bounce, splash

SourceFiles
    KinematicSurfaceFilm.C
    KinematicSurfaceFilmI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_KinematicSurfaceFilm_H
#define Foam_KinematicSurfaceFilm_H

#include "SurfaceFilmModel.H"
#include "UPtrList.H"
#include "liquidMixtureProperties.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
namespace regionModels
{
    namespace surfaceFilmModels
    {
        class surfaceFilmRegionModel;
    }
    namespace areaSurfaceFilmModels
    {
        class liquidFilmBase;
    }
}

/*---------------------------------------------------------------------------*\
                      Class KinematicSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class KinematicSurfaceFilm
:
    public SurfaceFilmModel<CloudType>
{
public:

    // Public Data

        //- Options for the interaction types
        enum class interactionType
        {
            absorb,       //!< absorb
            bounce,       //!< bounce
            splashBai     //!< Bai splash model
        };

        //- Names for interactionType
        static const Enum<interactionType> interactionTypeNames;


protected:

    // Data Types

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        typedef typename
            regionModels::areaSurfaceFilmModels::liquidFilmBase areaFilm;

        typedef typename
            regionModels::surfaceFilmModels::surfaceFilmRegionModel regionFilm;


    // Protected Data

        //- Reference to the cloud random number generator
        Random& rndGen_;


        // Region Film thermo

            //- Region Film liquid thermo
            liquidMixtureProperties* thermo_;

            //- Region Film reference pressure
            scalar pRef_;

            //- Region Film reference temperature
            scalar TRef_;

            //- Pointer to filmModel
            regionFilm* filmModel_;


        // Area Films

            //- UPointers to area films
            UPtrList<areaFilm> areaFilms_;


        // Interaction model data

            //- Interaction type enumeration
            interactionType interactionType_;

            //- Particle type IDs that can interact with the film
            //  If empty (default) all type IDs are included
            labelList parcelTypes_;

            //- Film thickness beyond which patch is assumed to be wet
            scalar deltaWet_;

            //- Splash parcel type label - id assigned to identify parcel for
            //  post-processing. If not specified, defaults to originating cloud
            //  type
            label splashParcelType_;

            //- Number of new parcels resulting from splash event
            label parcelsPerSplash_;

            //- Maximum splash particle diameter for Chi-square distribution
            //  Default is incident particle diameter
            scalar dMaxSplash_;

            //- Minimum splash particle diameter for Chi-square distribution
            //  Default is 0.1 dMaxSplash
            scalar dMinSplash_;

            // Surface roughness coefficient typically in the range 1300 - 5200
            // and decreases with increasing surface roughness

                //- Dry surface roughness coefficient
                //  = 2630 for dry interaction (ref. Bai)
                scalar Adry_;

                //- Wet surface roughness coefficient
                //  = 1320 for wet interaction (ref. Bai)
                scalar Awet_;


            //- Skin friction typically in the range 0.6 < Cf < 0.8
            scalar Cf_;

            //- Counter for number of new splash parcels
            label nParcelsSplashed_;


    // Protected Member Functions

        //- Return a vector tangential to input vector, v
        vector tangentVector(const vector& v) const;

        //- Return splashed parcel direction
        vector splashDirection
        (
            const vector& tanVec1,
            const vector& tanVec2,
            const vector& nf
        ) const;

        //- Initialise thermo
        void init(bool binitThermo);

        //- Initialise pointers of films
        void initFilmModels();


        // Injection from sheet (ejection) helper functions

            //- Cache the film fields in preparation for injection
            virtual void cacheFilmFields(const areaFilm& film);

            //- Cache the film fields in preparation for injection
            virtual void cacheFilmFields
            (
                const label filmPatchi,
                const label primaryPatchi,
                const regionModels::surfaceFilmModels::surfaceFilmRegionModel&
            );

            //- Set the individual parcel properties
            virtual void setParcelProperties
            (
                parcelType& p,
                const label filmFacei
            ) const;


public:

    //- Runtime type information
    TypeName("kinematicSurfaceFilm");


    // Constructors

        //- Construct from components
        KinematicSurfaceFilm
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type = typeName,
            bool initThermo = true
        );

        //- Construct copy
        KinematicSurfaceFilm
        (
            const KinematicSurfaceFilm<CloudType>& sfm,
            bool initThermo = true
        );

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<SurfaceFilmModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceFilmModel<CloudType>>
            (
                new KinematicSurfaceFilm<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~KinematicSurfaceFilm() = default;


    // Member Functions

        // Interaction models

            //- Absorb parcel into film
            template<class filmType>
            void absorbInteraction
            (
                filmType&,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                const scalar mass,
                bool& keepParticle
            );

            //- Bounce parcel (flip parcel normal velocity)
            void bounceInteraction
            (
                parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            ) const;

            //- Parcel interaction with dry surface
            template<class filmType>
            void drySplashInteraction
            (
                filmType&,
                const scalar sigma,
                const scalar mu,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            );

            //- Parcel interaction with wetted surface
            template<class filmType>
            void wetSplashInteraction
            (
                filmType&,
                const scalar sigma,
                const scalar mu,
                parcelType& p,
                const polyPatch& pp,
                const label facei,
                bool& keepParticle
            );

            //- Bai parcel splash interaction model
            template<class filmType>
            void splashInteraction
            (
                filmType&,
                const parcelType& p,
                const polyPatch& pp,
                const label facei,
                const scalar mRatio,
                const scalar We,
                const scalar Wec,
                const scalar sigma,
                bool& keepParticle
            );


        // Evaluation

            //- Transfer parcel from cloud to surface film
            //  Returns true if parcel is to be transferred
            virtual bool transferParcel
            (
                parcelType& p,
                const polyPatch& pp,
                bool& keepParticle
            );


        // I-O

            //- Write surface film info
            virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "KinematicSurfaceFilm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
