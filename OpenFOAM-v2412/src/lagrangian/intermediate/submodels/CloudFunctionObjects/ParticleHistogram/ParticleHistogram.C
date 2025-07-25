/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2023 OpenCFD Ltd.
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

#include "ParticleHistogram.H"
#include "Pstream.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleHistogram<CloudType>::writeFileHeader(Ostream& os) const
{
    this->writeHeaderValue(os, "nBin", nBins_);
    this->writeHeaderValue(os, "min", range_.min());
    this->writeHeaderValue(os, "max", range_.max());
    this->writeHeader(os, "");
    this->writeCommented(os, "dEdge1");
    os  << tab << "dEdge2"
        << tab << "nParticles"
        << tab << "nParticlesCumulative"
        << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleHistogram<CloudType>::ParticleHistogram
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    functionObjects::writeFile
    (
        owner,
        this->localPath(),
        typeName
    ),
    collector_(this->coeffDict(), owner.mesh()),
    nBins_
    (
        this->coeffDict().template getCheck<label>("nBins", labelMinMax::ge(1))
    ),
    maxStoredParcels_(this->coeffDict().getScalar("maxStoredParcels")),
    range_
    (
        this->coeffDict().getScalar("min"),
        this->coeffDict().getScalar("max")
    ),
    binEdges_(nBins_ + 1),
    nParticlesCumulative_(),
    dParticles_(),
    nParticles_()
{
    writeFile::read(this->coeffDict());

    if (!range_.good())
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "Invalid histogram range: " << range_
            << exit(FatalIOError);
    }

    if (maxStoredParcels_ <= 0)
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "maxStoredParcels = " << maxStoredParcels_
            << ", cannot be equal to or less than zero"
            << exit(FatalIOError);
    }

    // Compute histogram-bin properties
    binEdges_[0] = range_.min();
    const scalar delta = range_.span()/scalar(nBins_);
    for (label i = 0; i < nBins_; ++i)
    {
        const scalar next = range_.min() + (i+1)*delta;
        binEdges_[i+1] = next;
    }

    const label sz = collector_.size();
    nParticlesCumulative_ = List<scalarList>(sz, scalarList(nBins_, Zero));
    dParticles_.resize(sz);
    nParticles_.resize(sz);
}


template<class CloudType>
Foam::ParticleHistogram<CloudType>::ParticleHistogram
(
    const ParticleHistogram<CloudType>& ph
)
:
    CloudFunctionObject<CloudType>(ph),
    writeFile(ph),
    collector_(ph.collector_),
    nBins_(ph.nBins_),
    maxStoredParcels_(ph.maxStoredParcels_),
    range_(ph.range_),
    binEdges_(ph.binEdges_),
    nParticlesCumulative_(ph.nParticlesCumulative_),
    dParticles_(ph.dParticles_),
    nParticles_(ph.nParticles_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ParticleHistogram<CloudType>::postPatch
(
    const parcelType& p,
    const polyPatch& pp,
    const typename parcelType::trackingData& td
)
{
    if (!collector_.isPatch())
    {
        return true;
    }

    const label patchi = pp.index();
    const label localPatchi = collector_.IDs().find(patchi);

    if
    (
        localPatchi != -1
     && dParticles_[localPatchi].size() < maxStoredParcels_
    )
    {
        dParticles_[localPatchi].append(p.d());
        nParticles_[localPatchi].append(p.nParticle());
    }

    return true;
}


template<class CloudType>
bool Foam::ParticleHistogram<CloudType>::postFace
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    if (collector_.isPatch())
    {
        return true;
    }

    const labelList& IDs = collector_.IDs();
    const List<boundBox>& BBs = collector_.BBs();
    const faceZoneMesh& fzm = this->owner().mesh().faceZones();

    forAll(IDs, i)
    {
        if (!BBs[i].contains(p.position()))
        {
            // Quick reject if the particle is not in the face zone bound box
            continue;
        }

        const label zonei = IDs[i];
        const label localFacei = fzm[zonei].find(p.face());

        if
        (
            localFacei != -1
         && dParticles_[i].size() < maxStoredParcels_
        )
        {
            dParticles_[i].append(p.d());
            nParticles_[i].append(p.nParticle());
        }
    }

    return true;
}


template<class CloudType>
void Foam::ParticleHistogram<CloudType>::write()
{
    const wordList& names = collector_.names();

    forAll(names, i)
    {
        List<scalarList> procDiameters(Pstream::nProcs());
        procDiameters[Pstream::myProcNo()] = dParticles_[i];
        Pstream::gatherList(procDiameters);

        List<scalarList> procParticles(Pstream::nProcs());
        procParticles[Pstream::myProcNo()] = nParticles_[i];
        Pstream::gatherList(procParticles);

        if (Pstream::master())
        {
            scalarList globalDiameters;
            globalDiameters = ListListOps::combine<scalarList>
            (
                procDiameters,
                accessOp<scalarList>()
            );

            scalarList globalParticles;
            globalParticles = ListListOps::combine<scalarList>
            (
                procParticles,
                accessOp<scalarList>()
            );

            // Compute histogram
            scalarList nParticles(nBins_, Zero);
            const scalar delta = range_.span()/scalar(nBins_);
            forAll(globalDiameters, j)
            {
                const label bini = (globalDiameters[j] - range_.min())/delta;
                if (bini >= 0 && bini < nBins_)
                {
                    nParticles[bini] += globalParticles[j];
                    nParticlesCumulative_[i][bini] += globalParticles[j];
                }
            }

            if (this->writeToFile())
            {
                autoPtr<OFstream> osPtr = this->newFileAtTime
                (
                    names[i],
                    this->owner().time().value()
                );
                OFstream& os = osPtr.ref();

                writeFileHeader(os);

                forAll(nParticles, j)
                {
                    os  << binEdges_[j] << tab
                        << binEdges_[j + 1] << tab
                        << nParticles[j] << tab
                        << nParticlesCumulative_[i][j]
                        << nl;
                }
            }
        }

        dParticles_[i].clearStorage();
        nParticles_[i].clearStorage();
    }
}


// ************************************************************************* //
