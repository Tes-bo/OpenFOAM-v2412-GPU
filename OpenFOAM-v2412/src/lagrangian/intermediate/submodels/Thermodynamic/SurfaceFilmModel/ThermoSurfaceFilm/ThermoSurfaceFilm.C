/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "ThermoSurfaceFilm.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const dictionary& dict,
    CloudType& owner
)
:
    KinematicSurfaceFilm<CloudType>(dict, owner, typeName,  false),
    thermo_
    (
        owner.db().objectRegistry::template lookupObject<SLGThermo>("SLGThermo")
    ),
    TFilmPatch_(),
    CpFilmPatch_()
{}


template<class CloudType>
Foam::ThermoSurfaceFilm<CloudType>::ThermoSurfaceFilm
(
    const ThermoSurfaceFilm<CloudType>& sfm
)
:
    KinematicSurfaceFilm<CloudType>(sfm, false),
    thermo_(sfm.thermo_),
    TFilmPatch_(sfm.TFilmPatch_),
    CpFilmPatch_(sfm.CpFilmPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
template<class filmType>
void Foam::ThermoSurfaceFilm<CloudType>::absorbInteraction
(
    filmType& film,
    const parcelType& p,
    const polyPatch& pp,
    const label facei,
    const scalar mass,
    bool& keepParticle
)
{
    DebugInfo<< "Parcel " << p.origId() << " absorbInteraction" << endl;

    // Patch face normal
    const vector& nf = pp.faceNormals()[facei];

    // Patch velocity
    const vector& Up = this->owner().U().boundaryField()[pp.index()][facei];

    // Relative parcel velocity
    const vector Urel(p.U() - Up);

    // Parcel normal velocity
    const vector Un(nf*(Urel & nf));

    // Parcel tangential velocity
    const vector Ut(Urel - Un);

    film.addSources
    (
        pp.index(),
        facei,
        mass,                             // mass
        mass*Ut,                          // tangential momentum
        mass*mag(Un),                     // impingement pressure
        mass*p.hs()                       // energy
    );

    this->nParcelsTransferred()++;

    this->totalMassTransferred() += mass;

    keepParticle = false;
}


template<class CloudType>
bool Foam::ThermoSurfaceFilm<CloudType>::transferParcel
(
    parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    const label patchi = pp.index();
    const label meshFacei = p.face();
    const label facei = pp.whichFace(meshFacei);

    this->initFilmModels();

    // Check the singleLayer film models
    if (this->filmModel_ && this->filmModel_->isRegionPatch(patchi))
    {
        auto& film = *(this->filmModel_);

        switch (this->interactionType_)
        {
            case KinematicSurfaceFilm<CloudType>::interactionType::bounce:
            {
                this->bounceInteraction(p, pp, facei, keepParticle);

                break;
            }

            case KinematicSurfaceFilm<CloudType>::interactionType::absorb:
            {
                const scalar m = p.nParticle()*p.mass();

                this->absorbInteraction //<regionFilm>
                    (film, p, pp, facei, m, keepParticle);

                break;
            }

            case KinematicSurfaceFilm<CloudType>::interactionType::splashBai:
            {
                // Local pressure
                const scalar pc = thermo_.thermo().p()[p.cell()];
                const liquidProperties& liq = thermo_.liquids().properties()[0];
                const scalar sigma = liq.sigma(pc, p.T());
                const scalar mu = liq.mu(pc, p.T());

                const bool dry
                (
                    this->deltaFilmPatch_[patchi][facei] < this->deltaWet_
                );

                if (dry)
                {
                    this->drySplashInteraction //<CloudType, regionFilm>
                        (film, sigma, mu, p, pp, facei, keepParticle);
                }
                else
                {
                    this->wetSplashInteraction //<regionFilm>
                        (film, sigma, mu, p, pp, facei, keepParticle);
                }

                break;
            }

            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type enumeration"
                    << abort(FatalError);
            }
        }

        // Transfer parcel/parcel interactions complete
        return true;
    }


    // Check the area film models
    for (areaFilm& film : this->areaFilms_)
    {
        const label filmFacei
        (
            film.isRegionPatch(patchi)
          ? film.regionMesh().whichFace(meshFacei)
          : -1
        );

        if (filmFacei < 0)
        {
            // Film model does not include this patch face
            continue;
        }

        switch (this->interactionType_)
        {
            case KinematicSurfaceFilm<CloudType>::interactionType::bounce:
            {
                this->bounceInteraction(p, pp, facei, keepParticle);

                break;
            }

            case KinematicSurfaceFilm<CloudType>::interactionType::absorb:
            {
                const scalar m = p.nParticle()*p.mass();

                this->absorbInteraction //<areaFilm>
                (
                    film, p, pp, facei, m, keepParticle
                );
                break;
            }

            case KinematicSurfaceFilm<CloudType>::interactionType::splashBai:
            {
                // Local pressure
                const scalar pc = thermo_.thermo().p()[p.cell()];
                const liquidProperties& liq = thermo_.liquids().properties()[0];
                const scalar sigma = liq.sigma(pc, p.T());
                const scalar mu = liq.mu(pc, p.T());

                const bool dry
                (
                    film.h()[filmFacei] < this->deltaWet_
                );

                if (dry)
                {
                    this->drySplashInteraction //<areaFilm>
                        (film, sigma, mu, p, pp, facei, keepParticle);
                }
                else
                {
                    this->wetSplashInteraction //<areaFilm>
                        (film, sigma, mu, p, pp, facei, keepParticle);
                }

                break;
            }

            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type enumeration"
                    << abort(FatalError);
            }
        }

        // Transfer parcel/parcel interactions complete
        return true;
    }

    // Parcel did not interact with film
    return false;
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const label filmPatchi,
    const label primaryPatchi,
    const regionModels::surfaceFilmModels::surfaceFilmRegionModel& filmModel
)
{
    KinematicSurfaceFilm<CloudType>::cacheFilmFields
    (
        filmPatchi,
        primaryPatchi,
        filmModel
    );

    TFilmPatch_ = filmModel.Ts().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, TFilmPatch_);

    CpFilmPatch_ = filmModel.Cp().boundaryField()[filmPatchi];
    filmModel.toPrimary(filmPatchi, CpFilmPatch_);
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::cacheFilmFields
(
    const regionModels::areaSurfaceFilmModels::liquidFilmBase& film
)
{
    KinematicSurfaceFilm<CloudType>::cacheFilmFields(film);

    // Direct copy (one-to-one mapping)
    TFilmPatch_ = film.Tf().primitiveField();

    // Direct copy (one-to-one mapping)
    CpFilmPatch_ = film.Cp().primitiveField();
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::setParcelProperties
(
    parcelType& p,
    const label filmFacei
) const
{
    KinematicSurfaceFilm<CloudType>::setParcelProperties(p, filmFacei);

    // Set parcel properties
    p.T() = TFilmPatch_[filmFacei];
    p.Cp() = CpFilmPatch_[filmFacei];
}


template<class CloudType>
void Foam::ThermoSurfaceFilm<CloudType>::info()
{
    KinematicSurfaceFilm<CloudType>::info();
}


// ************************************************************************* //
