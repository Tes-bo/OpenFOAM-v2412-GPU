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

#include "LocalInteraction.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class CloudType>
void Foam::LocalInteraction<CloudType>::writeFileHeader(Ostream& os)
{
    PatchInteractionModel<CloudType>::writeFileHeader(os);

    forAll(nEscape_, patchi)
    {
        const word& patchName = patchData_[patchi].patchName();

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


template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    patchData_(cloud.mesh(), this->coeffDict()),
    nEscape_(patchData_.size()),
    massEscape_(nEscape_.size()),
    nStick_(nEscape_.size()),
    massStick_(nEscape_.size()),
    writeFields_(this->coeffDict().getOrDefault("writeFields", false)),
    injIdToIndex_(),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{
    const bool outputByInjectorId
        = this->coeffDict().getOrDefault("outputByInjectorId", false);

    if (writeFields_)
    {
        Info<< "    Interaction fields will be written to "
            << IOobject::scopedName(this->owner().name(), "massEscape")
            << " and "
            << IOobject::scopedName(this->owner().name(), "massStick") << endl;

        (void)massEscape();
        (void)massStick();
    }
    else
    {
        Info<< "    Interaction fields will not be written" << endl;
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
    if (!nInjectors)
    {
        nInjectors = 1;
    }

    // Check that interactions are valid/specified
    forAll(patchData_, patchi)
    {
        const word& interactionTypeName =
            patchData_[patchi].interactionTypeName();
        const typename PatchInteractionModel<CloudType>::interactionType& it =
            this->wordToInteractionType(interactionTypeName);

        if (it == PatchInteractionModel<CloudType>::itOther)
        {
            const word& patchName = patchData_[patchi].patchName();
            FatalErrorInFunction
                << "Unknown patch interaction type "
                << interactionTypeName << " for patch " << patchName
                << ". Valid selections are:"
                << this->PatchInteractionModel<CloudType>::interactionTypeNames_
                << nl << exit(FatalError);
        }

        nEscape_[patchi].setSize(nInjectors, Zero);
        massEscape_[patchi].setSize(nInjectors, Zero);
        nStick_[patchi].setSize(nInjectors, Zero);
        massStick_[patchi].setSize(nInjectors, Zero);
    }
}


template<class CloudType>
Foam::LocalInteraction<CloudType>::LocalInteraction
(
    const LocalInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    patchData_(pim.patchData_),
    nEscape_(pim.nEscape_),
    massEscape_(pim.massEscape_),
    nStick_(pim.nStick_),
    massStick_(pim.massStick_),
    writeFields_(pim.writeFields_),
    injIdToIndex_(pim.injIdToIndex_),
    massEscapePtr_(nullptr),
    massStickPtr_(nullptr)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massEscape()
{
    if (!massEscapePtr_)
    {
        const fvMesh& mesh = this->owner().mesh();

        massEscapePtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::scopedName(this->owner().name(), "massEscape"),
                    mesh.time().timeName(),
                    mesh.thisDb(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh,
                dimensionedScalar(dimMass, Zero)
            )
        );
    }

    return *massEscapePtr_;
}


template<class CloudType>
Foam::volScalarField& Foam::LocalInteraction<CloudType>::massStick()
{
    if (!massStickPtr_)
    {
        const fvMesh& mesh = this->owner().mesh();

        massStickPtr_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    IOobject::scopedName(this->owner().name(), "massStick"),
                    mesh.time().timeName(),
                    mesh.thisDb(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE,
                    IOobject::REGISTER
                ),
                mesh,
                dimensionedScalar(dimMass, Zero)
            )
        );
    }

    return *massStickPtr_;
}


template<class CloudType>
bool Foam::LocalInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    const label patchi = patchData_.applyToPatch(pp.index());

    if (patchi >= 0)
    {
        vector& U = p.U();

        // Location for storing the stats.
        const label idx =
        (
            injIdToIndex_.size()
          ? injIdToIndex_.lookup(p.typeId(), 0)
          : 0
        );

        typename PatchInteractionModel<CloudType>::interactionType it =
            this->wordToInteractionType
            (
                patchData_[patchi].interactionTypeName()
            );

        switch (it)
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

                const scalar dm = p.mass()*p.nParticle();

                nEscape_[patchi][idx]++;
                massEscape_[patchi][idx] += dm;

                if (writeFields_)
                {
                    const label pI = pp.index();
                    const label fI = pp.whichFace(p.face());
                    massEscape().boundaryFieldRef()[pI][fI] += dm;
                }
                break;
            }
            case PatchInteractionModel<CloudType>::itStick:
            {
                keepParticle = true;
                p.active(false);
                U = Zero;

                const scalar dm = p.mass()*p.nParticle();

                nStick_[patchi][idx]++;
                massStick_[patchi][idx] += dm;

                if (writeFields_)
                {
                    const label pI = pp.index();
                    const label fI = pp.whichFace(p.face());
                    massStick().boundaryFieldRef()[pI][fI] += dm;
                }
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
                    U -= (1.0 + patchData_[patchi].e())*Un*nw;
                }

                U -= patchData_[patchi].mu()*Ut;

                // Return velocity to global space
                U += Up;

                break;
            }
            default:
            {
                FatalErrorInFunction
                    << "Unknown interaction type "
                    << patchData_[patchi].interactionTypeName()
                    << "(" << it << ") for patch "
                    << patchData_[patchi].patchName()
                    << ". Valid selections are:" << this->interactionTypeNames_
                    << endl << abort(FatalError);
            }
        }

        return true;
    }

    return false;
}


template<class CloudType>
void Foam::LocalInteraction<CloudType>::info()
{
    PatchInteractionModel<CloudType>::info();

    // retrieve any stored data
    labelListList npe0(patchData_.size());
    scalarListList mpe0(patchData_.size());
    labelListList nps0(patchData_.size());
    scalarListList mps0(patchData_.size());

    forAll(patchData_, patchi)
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

    // accumulate current data
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
        // Since injIdToIndex_ is a one-to-one mapping (starting at zero),
        // can simply invert it.
        labelList indexToInjector(injIdToIndex_.size());
        forAllConstIters(injIdToIndex_, iter)
        {
            indexToInjector[iter.val()] = iter.key();
        }

        forAll(patchData_, patchi)
        {
            forAll(mpe[patchi], indexi)
            {
                const word& patchName = patchData_[patchi].patchName();

                Log_<< "    Parcel fate: patch " <<  patchName
                    << " (number, mass)" << nl
                    << "      - escape  (injector " << indexToInjector[indexi]
                    << " )  = " << npe[patchi][indexi]
                    << ", " << mpe[patchi][indexi] << nl
                    << "      - stick   (injector " << indexToInjector[indexi]
                    << " )  = " << nps[patchi][indexi]
                    << ", " << mps[patchi][indexi] << nl;
            }
        }
    }
    else
    {
        forAll(patchData_, patchi)
        {
            const word& patchName = patchData_[patchi].patchName();

            Log_<< "    Parcel fate: patch " << patchName
                << " (number, mass)" << nl
                << "      - escape                      = "
                << npe[patchi][0] << ", " << mpe[patchi][0] << nl
                << "      - stick                       = "
                << nps[patchi][0] << ", " << mps[patchi][0] << nl;
        }
    }

    forAll(npe, patchi)
    {
        forAll(npe[patchi], injectori)
        {
            this->file()
                << tab << npe[patchi][injectori]
                << tab << mpe[patchi][injectori]
                << tab << nps[patchi][injectori]
                << tab << mps[patchi][injectori];
        }
    }

    this->file() << endl;

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
