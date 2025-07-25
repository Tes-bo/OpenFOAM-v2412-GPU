/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "StandardWallInteraction.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::StandardWallInteraction<CloudType>::writeFileHeader(Ostream& os)
{
    PatchInteractionModel<CloudType>::writeFileHeader(os);

    forAll(nEscape_, patchi)
    {
        const word& patchName = mesh_.boundary()[patchi].name();

        forAll(nEscape_[patchi], injectori)
        {
            const word suffix = Foam::name(injectori);
            this->writeTabbed(os, patchName + "_nEscape_" + suffix);
            this->writeTabbed(os, patchName + "_massEscape_" + suffix);
            this->writeTabbed(os, patchName + "_nStick_" + suffix);
            this->writeTabbed(os, patchName + "_massStick_" + suffix);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    mesh_(cloud.mesh()),
    interactionType_
    (
        this->wordToInteractionType(this->coeffDict().getWord("type"))
    ),
    e_(0.0),
    mu_(0.0),
    nEscape_(mesh_.boundaryMesh().nNonProcessor()),
    massEscape_(nEscape_.size()),
    nStick_(nEscape_.size()),
    massStick_(nEscape_.size()),
    injIdToIndex_()
{
    const bool outputByInjectorId =
        this->coeffDict().getOrDefault("outputByInjectorId", false);

    switch (interactionType_)
    {
        case PatchInteractionModel<CloudType>::itOther:
        {
            const word interactionTypeName(this->coeffDict().getWord("type"));

            FatalErrorInFunction
                << "Unknown interaction result type "
                << interactionTypeName
                << ". Valid selections are:" << this->interactionTypeNames_
                << endl << exit(FatalError);

            break;
        }
        case PatchInteractionModel<CloudType>::itRebound:
        {
            e_ = this->coeffDict().getOrDefault("e", 1.0);
            mu_ = this->coeffDict().getOrDefault("mu", 0.0);
            break;
        }
        default:
        {}
    }

    // Determine the number of injectors and the injector mapping
    label nInjectors = 0;
    if (outputByInjectorId)
    {
        for (const auto& inj : cloud.injectors())
        {
            injIdToIndex_.insert(inj.injectorID(), nInjectors++);
        }
    }

    // The normal case, and safety if injector mapping was somehow null.
    if (injIdToIndex_.empty())
    {
        nInjectors = 1;
    }

    forAll(nEscape_, patchi)
    {
        nEscape_[patchi].setSize(nInjectors, Zero);
        massEscape_[patchi].setSize(nInjectors, Zero);
        nStick_[patchi].setSize(nInjectors, Zero);
        massStick_[patchi].setSize(nInjectors, Zero);
    }
}


template<class CloudType>
Foam::StandardWallInteraction<CloudType>::StandardWallInteraction
(
    const StandardWallInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    mesh_(pim.mesh_),
    interactionType_(pim.interactionType_),
    e_(pim.e_),
    mu_(pim.mu_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    injIdToIndex_(pim.injIdToIndex_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::StandardWallInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    vector& U = p.U();

    if (isA<wallPolyPatch>(pp))
    {
        // Location for storing the stats.
        const label idx =
        (
            injIdToIndex_.size()
          ? injIdToIndex_.lookup(p.typeId(), 0)
          : 0
        );

        switch (interactionType_)
        {
            case PatchInteractionModel<CloudType>::itNone:
            {
                return false;
            }
            case PatchInteractionModel<CloudType>::itEscape:
            {
                keepParticle = false;
                p.active(false);
                U = Zero;

                const scalar dm = p.nParticle()*p.mass();

                nEscape_[pp.index()][idx]++;
                massEscape_[pp.index()][idx] += dm;
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                p.active(false);
                U = Zero;

                const scalar dm = p.nParticle()*p.mass();

                nStick_[pp.index()][idx]++;
                massStick_[pp.index()][idx] += dm;
                break;
            }
            case PatchInteractionModel<CloudType>::itRebound:
            {
                keepParticle = true;
                p.active(true);

                vector nw;
                vector Up;

                this->owner().patchData(p, pp, nw, Up);

                // Calculate motion relative to patch velocity
                U -= Up;

                if (mag(Up) > 0 && mag(U) < this->Urmax())
                {
                    WarningInFunction
                        << "Particle U the same as patch "
                        << "    The particle has been removed" << nl << endl;

                    keepParticle = false;
                    p.active(false);
                    U = Zero;
                    break;
                }

                scalar Un = U & nw;
                vector Ut = U - Un*nw;

                if (Un > 0)
                {
                    U -= (1.0 + e_)*Un*nw;
                }

                U -= mu_*Ut;

                // Return velocity to global space
                U += Up;


                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type "
                    << this->interactionTypeToWord(interactionType_)
                    << "(" << interactionType_ << ")" << endl
                    << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::StandardWallInteraction<CloudType>::info()
{
    PatchInteractionModel<CloudType>::info();

    labelListList npe0(nEscape_.size());
    scalarListList mpe0(nEscape_.size());
    labelListList nps0(nEscape_.size());
    scalarListList mps0(nEscape_.size());

    forAll(nEscape_, patchi)
    {
        label lsd = nEscape_[patchi].size();
        npe0[patchi].setSize(lsd, Zero);
        mpe0[patchi].setSize(lsd, Zero);
        nps0[patchi].setSize(lsd, Zero);
        mps0[patchi].setSize(lsd, Zero);
    }

    this->getModelProperty("nEscape", npe0);
    this->getModelProperty("massEscape", mpe0);
    this->getModelProperty("nStick", nps0);
    this->getModelProperty("massStick", mps0);

    // Accumulate current data
    labelListList npe(nEscape_);

    forAll(npe, i)
    {
        Pstream::listCombineGather(npe[i], plusEqOp<label>());
        npe[i] = npe[i] + npe0[i];
    }

    scalarListList mpe(massEscape_);
    forAll(mpe, i)
    {
        Pstream::listCombineGather(mpe[i], plusEqOp<scalar>());
        mpe[i] = mpe[i] + mpe0[i];
    }

    labelListList nps(nStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(nps[i], plusEqOp<label>());
        nps[i] = nps[i] + nps0[i];
    }

    scalarListList mps(massStick_);
    forAll(nps, i)
    {
        Pstream::listCombineGather(mps[i], plusEqOp<scalar>());
        mps[i] = mps[i] + mps0[i];
    }

    if (injIdToIndex_.size())
    {
        // Since injIdToIndex_ is a one-to-one mapping (starting as zero),
        // can simply invert it.
        labelList indexToInjector(injIdToIndex_.size());
        forAllConstIters(injIdToIndex_, iter)
        {
            indexToInjector[iter.val()] = iter.key();
        }

        forAll(npe, patchi)
        {
            forAll(mpe[patchi], indexi)
            {
                const word& patchName = mesh_.boundary()[patchi].name() ;

                Log_<< "    Parcel fate: patch " <<  patchName
                    << " (number, mass)" << nl
                    << "      - escape  (injector " << indexToInjector[indexi]
                    << ")  = " << npe[patchi][indexi]
                    << ", " << mpe[patchi][indexi] << nl
                    << "      - stick   (injector " << indexToInjector[indexi]
                    << ")  = " << nps[patchi][indexi]
                    << ", " << mps[patchi][indexi] << nl;

                this->file()
                    << tab << npe[patchi][indexi] << tab << mpe[patchi][indexi]
                    << tab << nps[patchi][indexi] << tab << mps[patchi][indexi];
            }
        }

        this->file() << endl;
    }
    else
    {
        forAll(npe, patchi)
        {
            const word& patchName = mesh_.boundary()[patchi].name();

            Log_<< "    Parcel fate: patch (number, mass) "
                << patchName << nl
                << "      - escape                      = "
                << npe[patchi][0] << ", " << mpe[patchi][0] << nl
                << "      - stick                       = "
                << nps[patchi][0] << ", " << mps[patchi][0] << nl;

            this->file()
                << tab << npe[patchi][0] << tab << mpe[patchi][0]
                << tab << nps[patchi][0] << tab << mps[patchi][0];
        }

        this->file() << endl;
    }

    if (this->writeTime())
    {
        this->setModelProperty("nEscape", npe);
        this->setModelProperty("massEscape", mpe);
        this->setModelProperty("nStick", nps);
        this->setModelProperty("massStick", mps);

        nEscape_ = Zero;
        massEscape_ = Zero;
        nStick_ = Zero;
        massStick_ = Zero;
    }
}


// ************************************************************************* //
