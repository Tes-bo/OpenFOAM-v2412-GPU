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

#include "ParticlePostProcessing.H"
#include "Pstream.H"
#include "ListOps.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticlePostProcessing<CloudType>::writeFileHeader(Ostream& os) const
{
    this->writeCommented(os, "Time");
    os  << ' ' << "currentProc";

    if (!header_.empty())
    {
        os  << ' ' << header_;
    }

    os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticlePostProcessing<CloudType>::ParticlePostProcessing
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
    maxStoredParcels_(this->coeffDict().getScalar("maxStoredParcels")),
    header_(),
    fields_(),
    times_(),
    data_()
{
    writeFile::read(this->coeffDict());

    this->coeffDict().readIfPresent("fields", fields_);

    if (maxStoredParcels_ <= 0)
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "maxStoredParcels = " << maxStoredParcels_
            << ", cannot be equal to or less than zero"
            << exit(FatalIOError);
    }

    const label sz = collector_.size();
    times_.resize(sz);
    data_.resize(sz);
}


template<class CloudType>
Foam::ParticlePostProcessing<CloudType>::ParticlePostProcessing
(
    const ParticlePostProcessing<CloudType>& ppp
)
:
    CloudFunctionObject<CloudType>(ppp),
    writeFile(ppp),
    collector_(ppp.collector_),
    maxStoredParcels_(ppp.maxStoredParcels_),
    header_(ppp.header_),
    fields_(ppp.fields_),
    times_(ppp.times_),
    data_(ppp.data_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::ParticlePostProcessing<CloudType>::postPatch
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

    if (header_.empty())
    {
        OStringStream data;
        p.writeProperties(data, fields_, " ", true);
        header_ = data.str();
    }

    if (localPatchi != -1 && data_[localPatchi].size() < maxStoredParcels_)
    {
        times_[localPatchi].append(this->owner().time().value());

        OStringStream data;
        data<< Pstream::myProcNo();
        p.writeProperties(data, fields_, " ", false);

        data_[localPatchi].append(data.str());
    }

    return true;
}


template<class CloudType>
bool Foam::ParticlePostProcessing<CloudType>::postFace
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

    if (header_.empty())
    {
        OStringStream data;
        p.writeProperties(data, fields_, " ", true);
        header_ = data.str();
    }

    forAll(IDs, i)
    {
        if (!BBs[i].contains(p.position()))
        {
            // Quick reject if the particle is not in the face zone bound box
            continue;
        }

        const label zonei = IDs[i];
        const label localFacei = fzm[zonei].find(p.face());

        if (localFacei != -1 && data_[localFacei].size() < maxStoredParcels_)
        {
            times_[i].append(this->owner().time().value());

            OStringStream data;
            data<< Pstream::myProcNo();
            p.writeProperties(data, fields_, " ", false);

            data_[i].append(data.str());
        }
    }

    return true;
}


template<class CloudType>
void Foam::ParticlePostProcessing<CloudType>::write()
{
    const wordList& names = collector_.names();

    forAll(names, i)
    {
        List<scalarList> procTimes(Pstream::nProcs());
        procTimes[Pstream::myProcNo()] = times_[i];
        Pstream::gatherList(procTimes);

        List<List<string>> procData(Pstream::nProcs());
        procData[Pstream::myProcNo()] = data_[i];
        Pstream::gatherList(procData);

        Pstream::combineReduce
        (
            header_,
            [](string& x, const string& y)
            {
                if (y.size() > x.size())
                {
                    x = y;
                }
            }
        );

        if (Pstream::master())
        {
            List<string> globalData;
            globalData = ListListOps::combine<List<string>>
            (
                procData,
                accessOp<List<string>>()
            );

            scalarList globalTimes;
            globalTimes = ListListOps::combine<scalarList>
            (
                procTimes,
                accessOp<scalarList>()
            );

            if (this->writeToFile())
            {
                autoPtr<OFstream> osPtr = this->newFileAtTime
                (
                    names[i],
                    this->owner().time().value()
                );
                OFstream& os = osPtr.ref();

                writeFileHeader(os);

                const labelList indices(sortedOrder(globalTimes));
                forAll(globalTimes, j)
                {
                    const label datai = indices[j];

                    os  << globalTimes[datai] << tab
                        << globalData[datai].c_str()
                        << nl;
                }
            }
        }

        times_[i].clearStorage();
        data_[i].clearStorage();
    }
}


// ************************************************************************* //
