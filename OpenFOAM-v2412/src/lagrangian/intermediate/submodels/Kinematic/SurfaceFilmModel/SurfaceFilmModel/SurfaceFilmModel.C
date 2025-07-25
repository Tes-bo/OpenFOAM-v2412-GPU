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

\*---------------------------------------------------------------------------*/

#include "SurfaceFilmModel.H"
#include "mathematicalConstants.H"
#include "surfaceFilmRegionModel.H"
#include "liquidFilmBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel(CloudType& owner)
:
    CloudSubModelBase<CloudType>(owner),
    g_(owner.g()),
    ejectedParcelType_(0),
    injectionOffset_(1.1),
    minDiameter_(0),
    massParcelPatch_(),
    diameterParcelPatch_(),
    UFilmPatch_(),
    rhoFilmPatch_(),
    deltaFilmPatch_(),
    nParcelsTransferred_(0),
    nParcelsInjected_(0),
    totalMassTransferred_(0)
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const dictionary& dict,
    CloudType& owner,
    const word& type
)
:
    CloudSubModelBase<CloudType>(owner, dict, typeName, type),
    g_(owner.g()),
    ejectedParcelType_
    (
        this->coeffDict().template getOrDefault<label>("ejectedParcelType", -1)
    ),
    injectionOffset_
    (
        this->coeffDict().template getOrDefault<scalar>("injectionOffset", 1.1)
    ),
    minDiameter_
    (
        this->coeffDict().template getOrDefault<scalar>("minDiameter", -1)
    ),
    massParcelPatch_(),
    diameterParcelPatch_(),
    UFilmPatch_(),
    rhoFilmPatch_(),
    deltaFilmPatch_(owner.mesh().boundary().size()),
    nParcelsTransferred_(0),
    nParcelsInjected_(0),
    totalMassTransferred_()
{}


template<class CloudType>
Foam::SurfaceFilmModel<CloudType>::SurfaceFilmModel
(
    const SurfaceFilmModel<CloudType>& sfm
)
:
    CloudSubModelBase<CloudType>(sfm),
    g_(sfm.g_),
    ejectedParcelType_(sfm.ejectedParcelType_),
    injectionOffset_(sfm.injectionOffset_),
    minDiameter_(sfm.minDiameter_),
    massParcelPatch_(sfm.massParcelPatch_),
    diameterParcelPatch_(sfm.diameterParcelPatch_),
    UFilmPatch_(sfm.UFilmPatch_),
    rhoFilmPatch_(sfm.rhoFilmPatch_),
    deltaFilmPatch_(sfm.deltaFilmPatch_),
    nParcelsTransferred_(sfm.nParcelsTransferred_),
    nParcelsInjected_(sfm.nParcelsInjected_),
    totalMassTransferred_(sfm.totalMassTransferred_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class CloudTrackType>
void Foam::SurfaceFilmModel<CloudType>::injectParticles
(
    const label primaryPatchi,
    const labelUList& injectorCells,
    CloudTrackType& cloud
)
{
    const fvMesh& mesh = this->owner().mesh();

    const vectorField& Cf = mesh.C().boundaryField()[primaryPatchi];
    const vectorField& Sf = mesh.Sf().boundaryField()[primaryPatchi];

    // Cached sizes are identical to the patch size
    forAll(injectorCells, facei)
    {
        const label celli = injectorCells[facei];

        if (diameterParcelPatch_[facei] > 0)
        {
            const scalar offset =
            (
                injectionOffset_
              * max
                (
                    diameterParcelPatch_[facei],
                    deltaFilmPatch_[primaryPatchi][facei]
                )
            );

            const point pos = Cf[facei] - offset*normalised(Sf[facei]);

            // Create a new parcel
            parcelType* pPtr =
                new parcelType(this->owner().pMesh(), pos, celli);

            // Check/set new parcel thermo properties
            cloud.setParcelThermoProperties(*pPtr, 0.0);

            setParcelProperties(*pPtr, facei);

            if (pPtr->nParticle() > 0.001)
            {
                // Check new parcel properties
                cloud.checkParcelProperties(*pPtr, 0.0, false);

                // Add the new parcel to the cloud
                cloud.addParticle(pPtr);

                ++nParcelsInjected_;
            }
            else
            {
                // TODO: cache mass and re-distribute?
                delete pPtr;
            }
        }
    }
}


template<class CloudType>
template<class CloudTrackType>
void Foam::SurfaceFilmModel<CloudType>::injectParticles
(
    const UList<labelPair>& patchFaces,
    CloudTrackType& cloud
)
{
    const fvMesh& mesh = this->owner().mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    const auto& Cf = mesh.C().boundaryField();
    const auto& Sf = mesh.Sf().boundaryField();

    forAll(patchFaces, filmFacei)
    {
        const labelPair& patchAndFace = patchFaces[filmFacei];
        const label patchi = patchAndFace.first();
        const label facei = patchAndFace.second();

        if (patchi < 0) continue;  // extra safety

        const label celli = pbm[patchi].faceCells()[facei];

        if (diameterParcelPatch_[filmFacei] > 0)
        {
            const scalar offset =
                injectionOffset_ * max
                (
                    diameterParcelPatch_[filmFacei],
                    deltaFilmPatch_[patchi][facei]
                );

            const point pos =
            (
                Cf[patchAndFace]
              - offset * normalised(Sf[patchAndFace])
            );

            // Create a new parcel
            parcelType* pPtr =
                new parcelType(this->owner().pMesh(), pos, celli);

            // Check/set new parcel thermo properties
            cloud.setParcelThermoProperties(*pPtr, 0.0);

            setParcelProperties(*pPtr, filmFacei);

            if (pPtr->nParticle() > 0.001)
            {
                // Check new parcel properties
                cloud.checkParcelProperties(*pPtr, 0.0, false);

                // Add the new parcel to the cloud
                cloud.addParticle(pPtr);

                ++nParcelsInjected_;
            }
            else
            {
                // TODO: cache mass and re-distribute?
                delete pPtr;
            }
        }
    }
}


template<class CloudType>
template<class TrackCloudType>
void Foam::SurfaceFilmModel<CloudType>::inject(TrackCloudType& cloud)
{
    if (!this->active())
    {
        return;
    }

    const fvMesh& mesh = this->owner().mesh();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // Check the singleLayer type of films
    {
        const auto* filmPtr =
            mesh.time().objectRegistry::template findObject<regionFilm>
            (
                "surfaceFilmProperties"
            );

        if (filmPtr && filmPtr->active())
        {
            const auto& film = *filmPtr;
            const labelList& filmPatches = film.intCoupledPatchIDs();
            const labelList& primaryPatches = film.primaryPatchIDs();

            forAll(filmPatches, i)
            {
                const label filmPatchi = filmPatches[i];
                const label primaryPatchi = primaryPatches[i];

                cacheFilmFields(filmPatchi, primaryPatchi, film);

                injectParticles
                (
                    primaryPatchi,
                    pbm[primaryPatchi].faceCells(),  // injector cells
                    cloud
                );
            }
        }
    }

    // Check finite area films
    for
    (
        const areaFilm& regionFa
      : mesh.time().objectRegistry::template csorted<areaFilm>()
    )
    {
        if (regionFa.active())
        {
            auto& film = const_cast<areaFilm&>(regionFa);

            const List<labelPair>& patchFaces =
                film.regionMesh().whichPatchFaces();


            cacheFilmFields(film);

            injectParticles(patchFaces, cloud);

            forAll(patchFaces, filmFacei)
            {
                const label patchi = patchFaces[filmFacei].first();
                const label facei = patchFaces[filmFacei].second();

                if (diameterParcelPatch_[filmFacei] > 0)
                {
                    film.addSources
                    (
                        patchi,
                        facei,
                       -massParcelPatch_[filmFacei],// mass
                        Zero,                       // tangential momentum
                        Zero,                       // impingement
                        Zero                        // energy
                    );
                }
            }
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& filmModel
)
{
    massParcelPatch_ = filmModel.cloudMassTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, massParcelPatch_);

    diameterParcelPatch_ =
        filmModel.cloudDiameterTrans().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, diameterParcelPatch_, maxEqOp<scalar>());

    UFilmPatch_ = filmModel.Us().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, UFilmPatch_);

    rhoFilmPatch_ = filmModel.rho().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, rhoFilmPatch_);

    deltaFilmPatch_[primaryPatchi] =
        filmModel.delta().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, deltaFilmPatch_[primaryPatchi]);
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::cacheFilmFields
(
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film
)
{
    const polyBoundaryMesh& pbm = this->owner().mesh().boundaryMesh();
    const volSurfaceMapping& map = film.region().vsm();

    // The polyPatch/local-face for each of the faceLabels
    const List<labelPair>& patchFaces =
        film.regionMesh().whichPatchFaces();

    const label nFaces = film.Uf().size();  // or regionMesh().nFaces()


    // Flat fields

    massParcelPatch_.resize(nFaces, Zero);
    map.mapToSurface(film.cloudMassTrans(), massParcelPatch_);

    diameterParcelPatch_.resize(nFaces, Zero);
    map.mapToSurface(film.cloudDiameterTrans(), diameterParcelPatch_);

    // Direct copy (one-to-one mapping)
    UFilmPatch_ = film.Uf().primitiveField();

    // Direct copy (one-to-one mapping)
    rhoFilmPatch_ = film.rho().primitiveField();


    // Per-patch fields

    // Same as film.region().primaryPatchIDs()
    for (const label patchi : film.regionMesh().whichPolyPatches())
    {
        deltaFilmPatch_[patchi].resize(pbm[patchi].size(), Zero);
    }

    forAll(patchFaces, i)
    {
        const label patchi = patchFaces[i].first();
        const label facei = patchFaces[i].second();

        if (patchi >= 0)
        {
            deltaFilmPatch_[patchi][facei] = film.h()[i];
        }
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    // Set parcel properties
    scalar vol = mathematical::pi/6.0*pow3(diameterParcelPatch_[filmFacei]);
    p.d() = diameterParcelPatch_[filmFacei];
    p.U() = UFilmPatch_[filmFacei];
    p.rho() = rhoFilmPatch_[filmFacei];

    p.nParticle() = massParcelPatch_[filmFacei]/p.rho()/vol;

    if (minDiameter_ != -1)
    {
        if (p.d() < minDiameter_)
        {
            p.nParticle() = 0;
        }
    }

    if (ejectedParcelType_ >= 0)
    {
        p.typeId() = ejectedParcelType_;
    }
}


template<class CloudType>
void Foam::SurfaceFilmModel<CloudType>::info()
{
    CloudSubModelBase<CloudType>::info();

    label nTrans0 =
        this->template getModelProperty<label>("nParcelsTransferred");

    label nInject0 =
        this->template getModelProperty<label>("nParcelsInjected");

    scalar massTransferred0 =
        this->template getModelProperty<scalar>("massTransferred");

    label nTransTotal =
        nTrans0 + returnReduce(nParcelsTransferred_, sumOp<label>());

    label nInjectTotal =
        nInject0 + returnReduce(nParcelsInjected_, sumOp<label>());

    scalar massTransferredTotal =
        massTransferred0 + returnReduce(totalMassTransferred_, sumOp<scalar>());


    Log_<< "    Surface film:" << nl
        << "      - parcels absorbed            = " << nTransTotal << nl
        << "      - mass absorbed               = " << massTransferredTotal << nl
        << "      - parcels ejected             = " << nInjectTotal << endl;

    if (this->writeTime())
    {
        this->setModelProperty("nParcelsTransferred", nTransTotal);
        this->setModelProperty("nParcelsInjected", nInjectTotal);
        this->setModelProperty("massTransferred", massTransferredTotal);

        nParcelsTransferred_ = 0;
        nParcelsInjected_ = 0;
        totalMassTransferred_ = 0;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SurfaceFilmModelNew.C"

// ************************************************************************* //
