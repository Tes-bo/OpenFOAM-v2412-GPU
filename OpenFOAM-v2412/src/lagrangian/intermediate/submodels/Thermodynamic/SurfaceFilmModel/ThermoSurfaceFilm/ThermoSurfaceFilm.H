/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021-2023 OpenCFD Ltd.
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
    Foam::ThermoSurfaceFilm

Group
    grpLagrangianIntermediateSurfaceFilmSubModels

Description
    Thermo parcel surface film model.

SourceFiles
    ThermoSurfaceFilm.C
    ThermoSurfaceFilmI.H

\*---------------------------------------------------------------------------*/

#ifndef Foam_ThermoSurfaceFilm_H
#define Foam_ThermoSurfaceFilm_H

#include "KinematicSurfaceFilm.H"
#include "UPtrList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class ThermoSurfaceFilm Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ThermoSurfaceFilm
:
    public KinematicSurfaceFilm<CloudType>
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

        //- Reference to the cloud thermo package
        const SLGThermo& thermo_;


        // Cached injector fields per film patch

            //- Film temperature / patch face
            scalarField TFilmPatch_;

            //- Film specific heat capacity / patch face
            scalarField CpFilmPatch_;


    // Protected Member Functions

        // Injection from sheet (ejection) helper functions

            //- Cache the film fields in preparation for injection
            virtual void cacheFilmFields
            (
                const label filmPatchi,
                const label primaryPatchi,
                const regionFilm&
            );

            //- Cache the finite area film fields in preparation for injection
            virtual void cacheFilmFields(const areaFilm&);

            //- Set the individual parcel properties
            virtual void setParcelProperties
            (
                parcelType& p,
                const label filmFacei
            ) const;


public:

    //- Runtime type information
    TypeName("thermoSurfaceFilm");


    // Constructors

        //- Construct from components
        ThermoSurfaceFilm(const dictionary& dict, CloudType& owner);

        //- Construct copy
        ThermoSurfaceFilm(const ThermoSurfaceFilm<CloudType>& sfm);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<SurfaceFilmModel<CloudType>> clone() const
        {
            return autoPtr<SurfaceFilmModel<CloudType>>
            (
                new ThermoSurfaceFilm<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~ThermoSurfaceFilm() = default;


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
    #include "ThermoSurfaceFilm.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
