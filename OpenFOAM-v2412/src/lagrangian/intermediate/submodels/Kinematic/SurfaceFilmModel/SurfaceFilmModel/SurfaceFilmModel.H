/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020-2022 OpenCFD Ltd.
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
    Foam::SurfaceFilmModel

Group
    grpLagrangianIntermediateSurfaceFilmSubModels

Description
    Templated wall surface film model class.

SourceFiles
    SurfaceFilmModel.C
    SurfaceFilmModelNew.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_SurfaceFilmModel_H
#define Foam_SurfaceFilmModel_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"


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
                      Class SurfaceFilmModel Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class SurfaceFilmModel
:
    public CloudSubModelBase<CloudType>
{
protected:

    // Data Types

        //- Convenience typedef to the cloud's parcel type
        typedef typename CloudType::parcelType parcelType;

        typedef typename
            regionModels::areaSurfaceFilmModels::liquidFilmBase
            areaFilm;

        typedef typename
            regionModels::surfaceFilmModels::surfaceFilmRegionModel
            regionFilm;


    // Protected Data

        //- Gravitational acceleration constant
        const dimensionedVector& g_;

        //- Ejected parcel type label - id assigned to identify parcel for
        //- post-processing. If not specified, defaults to originating cloud
        //- type
        label ejectedParcelType_;

        //- Injection offset position
        scalar injectionOffset_;

        //- Minimum diameter particle injection
        scalar minDiameter_;


        // Injector fields transferred from film patch

            //- Parcel mass / patch face
            scalarField massParcelPatch_;

            //- Parcel diameter / patch face
            scalarField diameterParcelPatch_;

            //- Film velocity / patch face
            Field<vector> UFilmPatch_;

            //- Film density / patch face
            scalarField rhoFilmPatch_;

            //- Film height of all film patches / patch face
            Field<scalarField> deltaFilmPatch_;


        // Counters

            //- Number of parcels transferred to the film model
            label nParcelsTransferred_;

            //- Number of parcels injected from the film model
            label nParcelsInjected_;


        // Total mass info

            //- Total mass transferred to the film
            scalar totalMassTransferred_;


    // Protected Functions

        //- Cache the film fields in preparation for injection
        virtual void cacheFilmFields
        (
            const label filmPatchi,
            const label primaryPatchi,
            const regionFilm&
        );

        //- Cache the finite area film fields in preparation for injection
        virtual void cacheFilmFields(const areaFilm& film);

        //- Inject particles in cloud
        template<class TrackCloudType>
        void injectParticles
        (
            const label primaryPatchi,
            const labelUList& injectorCells,  // patch-based
            TrackCloudType& cloud
        );

        //- Inject particles in cloud
        template<class TrackCloudType>
        void injectParticles
        (
            const UList<labelPair>& patchFaces,
            TrackCloudType& cloud
        );

        //- Set the individual parcel properties
        virtual void setParcelProperties
        (
            parcelType& p,
            const label filmFacei
        ) const;


public:

    //- Runtime type information
    TypeName("surfaceFilmModel");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        SurfaceFilmModel,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner
        ),
        (dict, owner)
    );


    // Constructors

        //- Construct null from owner
        SurfaceFilmModel(CloudType& owner);

        //- Construct from components
        SurfaceFilmModel
        (
            const dictionary& dict,
            CloudType& owner,
            const word& type
        );

        //- Construct copy
        SurfaceFilmModel(const SurfaceFilmModel<CloudType>& sfm);

        //- Construct and return a clone
        virtual autoPtr<SurfaceFilmModel<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~SurfaceFilmModel() = default;


    //- Selector
    static autoPtr<SurfaceFilmModel<CloudType>> New
    (
        const dictionary& dict,
        CloudType& owner
    );


    // Member Functions

    // Access

        //- Return gravitational acceleration constant
        inline const dimensionedVector& g() const noexcept;

        //- The number of parcels transferred to the film model
        inline label nParcelsTransferred() const noexcept;

        //- Non-const access to number of parcels transferred to the film model
        inline label& nParcelsTransferred() noexcept;

        //- The number of parcels injected from the film model
        inline label nParcelsInjected() const noexcept;

        //- Non-const access to number of parcels injected from the film model
        inline label& nParcelsInjected() noexcept;

        //- The total mass transferred
        inline scalar totalMassTransferred() const noexcept;

        //- Non-const access to the total mass transferred
        inline scalar& totalMassTransferred() noexcept;


    // Member Functions

        //- Transfer parcel from cloud to surface film
        //  Returns true if parcel is to be transferred
        virtual bool transferParcel
        (
            parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        ) = 0;

        //- Inject parcels into the cloud
        template<class TrackCloudType>
        void inject(TrackCloudType& cloud);


    // I-O

        //- Write surface film info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeSurfaceFilmModel(CloudType)                                        \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug                                        \
    (                                                                          \
        Foam::SurfaceFilmModel<kinematicCloudType>,                            \
        0                                                                      \
    );                                                                         \
    namespace Foam                                                             \
    {                                                                          \
        defineTemplateRunTimeSelectionTable                                    \
        (                                                                      \
            SurfaceFilmModel<kinematicCloudType>,                              \
            dictionary                                                         \
        );                                                                     \
    }


#define makeSurfaceFilmModelType(SS, CloudType)                                \
                                                                               \
    typedef Foam::CloudType::kinematicCloudType kinematicCloudType;            \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<kinematicCloudType>, 0);      \
                                                                               \
    Foam::SurfaceFilmModel<kinematicCloudType>::                               \
        adddictionaryConstructorToTable<Foam::SS<kinematicCloudType>>          \
            add##SS##CloudType##kinematicCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "SurfaceFilmModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
