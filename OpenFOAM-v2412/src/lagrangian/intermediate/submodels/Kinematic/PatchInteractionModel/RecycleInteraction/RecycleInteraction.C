/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

#include "RecycleInteraction.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::RecycleInteraction<CloudType>::writeFileHeader(Ostream& os)
{
    PatchInteractionModel<CloudType>::writeFileHeader(os);

    forAll(nRemoved_, i)
    {
        const word& outPatchName = recyclePatches_[i].first();

        forAll(nRemoved_[i], injectori)
        {
            const word suffix = Foam::name(injectori);
            this->writeTabbed(os, outPatchName + "_nRemoved_" + suffix);
            this->writeTabbed(os, outPatchName + "_massRemoved_" + suffix);
        }

        const word& inPatchName = recyclePatches_[i].second();

        forAll(nInjected_[i], injectori)
        {
            const word suffix = Foam::name(injectori);
            this->writeTabbed(os, inPatchName + "_nInjected_" + suffix);
            this->writeTabbed(os, inPatchName + "_massInjected_" + suffix);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RecycleInteraction<CloudType>::RecycleInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    mesh_(cloud.mesh()),
    recyclePatches_(this->coeffDict().lookup("recyclePatches")),
    recyclePatchesIds_(recyclePatches_.size()),
    recycledParcels_(recyclePatches_.size()),
    nRemoved_(recyclePatches_.size()), // per patch the no. of parcels
    massRemoved_(nRemoved_.size()),
    nInjected_(nRemoved_.size()),
    massInjected_(nRemoved_.size()),
    injectionPatchPtr_(nRemoved_.size()),
    recycleFraction_
    (
        this->coeffDict().template getCheck<scalar>
        (
            "recycleFraction",
            scalarMinMax::zero_one()
        )
    ),
    outputByInjectorId_
    (
        this->coeffDict().getOrDefault("outputByInjectorId", false)
    )
{
    // Determine the number of injectors and the injector mapping
    label nInjectors = 0;
    if (outputByInjectorId_)
    {
        for (const auto& inj : cloud.injectors())
        {
            injIdToIndex_.insert(inj.injectorID(), ++nInjectors);
        }
    }

    // The normal case, and safety if injector mapping was somehow null
    if (injIdToIndex_.empty())
    {
        nInjectors = 1;
    }

    forAll(nRemoved_, i)
    {
        // Create injection helper for each inflow patch
        injectionPatchPtr_.set
        (
            i,
            new patchInjectionBase(mesh_, recyclePatches_[i].second())
        );

        // Mappings from patch names to patch IDs
        recyclePatchesIds_[i].first() =
            mesh_.boundaryMesh().findPatchID(recyclePatches_[i].first());
        recyclePatchesIds_[i].second() =
            mesh_.boundaryMesh().findPatchID(recyclePatches_[i].second());

        // Storage for reporting
        nRemoved_[i].setSize(nInjectors, Zero);
        massRemoved_[i].setSize(nInjectors, Zero);
        nInjected_[i].setSize(nInjectors, Zero);
        massInjected_[i].setSize(nInjectors, Zero);
    }
}


template<class CloudType>
Foam::RecycleInteraction<CloudType>::RecycleInteraction
(
    const RecycleInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    mesh_(pim.mesh_),
    recyclePatches_(pim.recyclePatches_),
    recyclePatchesIds_(pim.recyclePatchesIds_),
    nRemoved_(pim.nRemoved_),
    massRemoved_(pim.massRemoved_),
    nInjected_(pim.nInjected_),
    massInjected_(pim.massInjected_),
    injIdToIndex_(pim.injIdToIndex_),
    injectionPatchPtr_(),
    recycleFraction_(pim.recycleFraction_),
    outputByInjectorId_(pim.outputByInjectorId_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::RecycleInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Injector ID
    const label idx =
    (
        injIdToIndex_.size()
      ? injIdToIndex_.lookup(p.typeId(), 0)
      : 0
    );

    // See if this patch is designated an outflow patch
    label addri = -1;
    forAll(recyclePatchesIds_, i)
    {
        if (recyclePatchesIds_[i].first() == pp.index())
        {
            addri = i;
            break;
        }
    }

    if (addri == -1)
    {
        // Not processing this outflow patch
        keepParticle = true;
        return false;
    }

    // Flag to remove current parcel and copy to local storage
    keepParticle = false;
    recycledParcels_[addri].append
    (
        static_cast<parcelType*>(p.clone().ptr())
    );

    ++nRemoved_[addri][idx];
    massRemoved_[addri][idx] += p.nParticle()*p.mass();

    return true;
}


template<class CloudType>
void Foam::RecycleInteraction<CloudType>::postEvolve()
{
    if (Pstream::parRun())
    {
        // See comments in Cloud::move() about transfer particles handling

        // Allocate transfer buffers
        PstreamBuffers pBufs;

        // Cache of opened UOPstream wrappers
        PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

        auto& rnd = this->owner().rndGen();

        forAll(recycledParcels_, addri)
        {
            auto& patchParcels = recycledParcels_[addri];
            auto& injectionPatch = injectionPatchPtr_[addri];

            for (parcelType& p : patchParcels)
            {
                // Choose a random location to insert the parcel
                const scalar fraction01 = rnd.template sample01<scalar>();

                // Identify the processor that owns the location
                const label toProci = injectionPatch.whichProc(fraction01);

                // Get/create output stream
                auto* osptr = UOPstreamPtrs.get(toProci);
                if (!osptr)
                {
                    osptr = new UOPstream(toProci, pBufs);
                    UOPstreamPtrs.set(toProci, osptr);
                }

                // Tuple: (address fraction particle)
                (*osptr) << addri << fraction01 << p;

                // Can now remove from list and delete
                delete(patchParcels.remove(&p));
            }
        }

        pBufs.finishedSends();

        // Not looping, so no early exit needed
        //
        // if (!returnReduceOr(pBufs.hasRecvData()))
        // {
        //     // No parcels to recycle
        //     return;
        // }

        // Retrieve from receive buffers
        for (const int proci : pBufs.allProcs())
        {
            if (pBufs.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs);

                // Read out each (address fraction particle) tuple
                while (!is.eof())
                {
                    const label addri = pTraits<label>(is);
                    const scalar fraction01 = pTraits<scalar>(is);
                    auto* newp = new parcelType(this->owner().mesh(), is);

                    // Parcel to be recycled
                    vector newPosition;
                    label cellOwner;
                    label dummy;
                    injectionPatchPtr_[addri].setPositionAndCell
                    (
                        mesh_,
                        fraction01,
                        this->owner().rndGen(),
                        newPosition,
                        cellOwner,
                        dummy,
                        dummy
                    );
                    newp->relocate(newPosition, cellOwner);
                    newp->nParticle() *= recycleFraction_;

                    // Assume parcel velocity is same as the carrier velocity
                    newp->U() = this->owner().U()[cellOwner];

                    // Injector ID
                    const label idx =
                    (
                        injIdToIndex_.size()
                      ? injIdToIndex_.lookup(newp->typeId(), 0)
                      : 0
                    );

                    ++nInjected_[addri][idx];
                    massInjected_[addri][idx] += newp->nParticle()*newp->mass();

                    this->owner().addParticle(newp);
                }
            }
        }
    }
    else
    {
        forAll(recycledParcels_, addri)
        {
            for (parcelType& p : recycledParcels_[addri])
            {
                parcelType* newp = recycledParcels_[addri].remove(&p);

                // Parcel to be recycled
                vector newPosition;
                label cellOwner;
                label dummy;
                injectionPatchPtr_[addri].setPositionAndCell
                (
                    mesh_,
                    this->owner().rndGen(),
                    newPosition,
                    cellOwner,
                    dummy,
                    dummy
                );

                // Update parcel properties
                newp->relocate(newPosition, cellOwner);
                newp->nParticle() *= recycleFraction_;

                // Assume parcel velocity is same as the carrier velocity
                newp->U() = this->owner().U()[cellOwner];

                // Injector ID
                const label idx =
                (
                    injIdToIndex_.size()
                  ? injIdToIndex_.lookup(newp->typeId(), 0)
                  : 0
                );
                ++nInjected_[addri][idx];
                massInjected_[addri][idx] += newp->nParticle()*newp->mass();

                this->owner().addParticle(newp);
            }
        }
    }
}


template<class CloudType>
void Foam::RecycleInteraction<CloudType>::info()
{
    PatchInteractionModel<CloudType>::info();

    labelListList npr0(nRemoved_.size());
    scalarListList mpr0(massRemoved_.size());
    labelListList npi0(nInjected_.size());
    scalarListList mpi0(massInjected_.size());

    forAll(nRemoved_, patchi)
    {
        const label lsd = nRemoved_[patchi].size();
        npr0[patchi].setSize(lsd, Zero);
        mpr0[patchi].setSize(lsd, Zero);
        npi0[patchi].setSize(lsd, Zero);
        mpi0[patchi].setSize(lsd, Zero);
    }

    this->getModelProperty("nRemoved", npr0);
    this->getModelProperty("massRemoved", mpr0);
    this->getModelProperty("nInjected", npi0);
    this->getModelProperty("massInjected", mpi0);

    // Accumulate current data
    labelListList npr(nRemoved_);

    forAll(npr, i)
    {
        Pstream::listCombineGather(npr[i], plusEqOp<label>());
        npr[i] = npr[i] + npr0[i];
    }

    scalarListList mpr(massRemoved_);
    forAll(mpr, i)
    {
        Pstream::listCombineGather(mpr[i], plusEqOp<scalar>());
        mpr[i] = mpr[i] + mpr0[i];
    }

    labelListList npi(nInjected_);
    forAll(npi, i)
    {
        Pstream::listCombineGather(npi[i], plusEqOp<label>());
        npi[i] = npi[i] + npi0[i];
    }

    scalarListList mpi(massInjected_);
    forAll(mpi, i)
    {
        Pstream::listCombineGather(mpi[i], plusEqOp<scalar>());
        mpi[i] = mpi[i] + mpi0[i];
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

        forAll(npr, i)
        {
            const word& outPatchName =  recyclePatches_[i].first();

            Log_<< "    Parcel fate: patch " <<  outPatchName
                << " (number, mass)" << nl;

            forAll(mpr[i], indexi)
            {
                Log_<< "      - removed  (injector " << indexToInjector[indexi]
                    << ")  = " << npr[i][indexi]
                    << ", " << mpr[i][indexi] << nl;

                this->file()
                    << tab << npr[i][indexi] << tab << mpr[i][indexi];
            }

            const word& inPatchName = recyclePatches_[i].second();

            Log_<< "    Parcel fate: patch " <<  inPatchName
                << " (number, mass)" << nl;

            forAll(mpi[i], indexi)
            {
                Log_<< "      - injected  (injector " << indexToInjector[indexi]
                    << ")  = " << npi[i][indexi]
                    << ", " << mpi[i][indexi] << nl;
                this->file()
                    << tab << npi[i][indexi] << tab << mpi[i][indexi];
            }
        }

        this->file() << endl;
    }
    else
    {
        forAll(npr, i)
        {
            const word& outPatchName =  recyclePatches_[i].first();

            Log_<< "    Parcel fate: patch " <<  outPatchName
                << " (number, mass)" << nl
                << "      - removed    = " << npr[i][0] << ", " << mpr[i][0]
                << nl;

            this->file()
                << tab << npr[i][0] << tab << mpr[i][0];
        }

        forAll(npi, i)
        {
            const word& inPatchName = recyclePatches_[i].second();

            Log_<< "    Parcel fate: patch " <<  inPatchName
                << " (number, mass)" << nl
                << "      - injected   = " << npi[i][0] << ", " << mpi[i][0]
                << nl;

            this->file()
                << tab << npi[i][0] << tab << mpi[i][0];
        }

        this->file() << endl;
    }

    if (this->writeTime())
    {
        this->setModelProperty("nRemoved", npr);
        this->setModelProperty("massRemoved", mpr);
        this->setModelProperty("nInjected", npi);
        this->setModelProperty("massInjected", mpi);

        nRemoved_ = Zero;
        massRemoved_ = Zero;
        nInjected_ = Zero;
        massInjected_ = Zero;
    }
}


// ************************************************************************* //
