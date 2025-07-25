/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 OpenCFD Ltd.
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

#include "ParticleZoneInfo.H"
#include "DynamicField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct particleInfoCombineOp
{
    void operator()(particleInfo& p1, const particleInfo& p2) const
    {
        // p2 not set
        if (p2.origID == -1)
        {
            return;
        }

        // p1 not set - initialise with p2
        if (p1.origID == -1)
        {
            p1 = p2;
            return;
        }

        // Set initial values
        if (p2.time0 < p1.time0)
        {
            p1.time0 = p2.time0;
            p1.d0 = p2.d0;
            p1.mass0 = p2.mass0;
        }

        // Accumulate age
        p1.age += p2.age;

        // Set latest available values
        if (p2.isOlderThan(p1))
        {
            p1.position = p2.position;
            p1.d = p2.d;
            p1.mass = p2.mass;
        }
    }
};

template<class Type>
Field<Type> getData
(
    const Foam::UList<Foam::particleInfo>& data,
    Type Foam::particleInfo::* field
)
{
    Field<Type> result(data.size());

    forAll(data, i)
    {
        result[i] = data[i].*field;
    }

    return result;
}

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleZoneInfo<CloudType>::writeWriter
(
    const UList<particleInfo>& data
)
{
    coordSet coords
    (
        "zoneParticles",
        "xyz",
        getData(data, &particleInfo::position),
        scalarList(data.size(), Zero)
    );

    writerPtr_->open(coords, this->baseTimeDir() / "zoneParticles");
    writerPtr_->beginTime(this->owner().time());

#undef  writeLocal
#define writeLocal(field)                                                      \
    writerPtr_->write(#field, getData(data, &particleInfo::field));

    writeLocal(origID);
    writeLocal(origProc);
    writeLocal(time0);
    writeLocal(age);
    writeLocal(d0);
    writeLocal(d);
    writeLocal(mass0);
    writeLocal(mass);
#undef  writeLocal

    writerPtr_->endTime();
    writerPtr_->close();
}


template<class CloudType>
void Foam::ParticleZoneInfo<CloudType>::writeFileHeader(Ostream& os) const
{
    this->writeHeaderValue(os, "cellZone", cellZoneName_);
    this->writeHeaderValue(os, "time", this->owner().time().timeOutputValue());
    this->writeHeader(os, "");
    this->writeCommented(os, "origID");
    os  << tab << "origProc"
        << tab << "(x y z)"
        << tab << "time0"
        << tab << "age"
        << tab << "d0"
        << tab << "d"
        << tab << "mass0"
        << tab << "mass"
        << endl;
}


template<class CloudType>
bool Foam::ParticleZoneInfo<CloudType>::inZone(const label celli) const
{
    return this->owner().mesh().cellZones()[cellZoneId_].whichCell(celli) != -1;
}


template<class CloudType>
Foam::label Foam::ParticleZoneInfo<CloudType>::getParticleID
(
    const particleInfo& p
) const
{
    forAll(data_, i)
    {
        const auto& d = data_[i];
        if ((d.origProc == p.origProc) && (d.origID == p.origID))
        {
            return i;
        }
    }

    return -1;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleZoneInfo<CloudType>::ParticleZoneInfo
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
    cellZoneName_(this->coeffDict().getWord("cellZone")),
    cellZoneId_(-1),
    data_(),
    movedParticles_(),
    maxIDs_(Pstream::nProcs(), Zero),
    writerPtr_(nullptr)
{
    if (Pstream::master())
    {
        const word writerType = this->coeffDict().getWord("writer");

        writerPtr_ = coordSetWriter::New
        (
            writerType,
            coordSetWriter::formatOptions(this->coeffDict(), writerType)
        );
    }

    writeFile::read(this->coeffDict());

    const auto& cellZones = owner.mesh().cellZones();

    cellZoneId_ = cellZones.findZoneID(cellZoneName_);
    if (cellZoneId_ == -1)
    {
        FatalIOErrorInFunction(this->coeffDict())
            << "Unable to find cellZone " << cellZoneName_
            << ". Available cellZones are:" << cellZones.names()
            << exit(FatalIOError);
    }

    Info<< "        Processing cellZone" << cellZoneName_ << " with id "
        << cellZoneId_ << endl;

    if (Pstream::master())
    {
        // Data was reduced on write
        labelList maxIDs;

        if (this->getModelProperty("maxIDs", maxIDs))
        {
            if (maxIDs.size() == Pstream::nProcs())
            {
                maxIDs_ = maxIDs;
                this->getModelProperty("data", data_);

                Info<< "        Restarting with " << data_.size()
                    << " particles" << endl;
            }
            else
            {
                WarningInFunction
                    << "Case restarted with a different number of processors."
                    << " Restarting particle statistics." << endl;

                // TODO
                // - use Cloud for base storage instead of local particle list?
            }
        }
    }
}


template<class CloudType>
Foam::ParticleZoneInfo<CloudType>::ParticleZoneInfo
(
    const ParticleZoneInfo<CloudType>& pzi
)
:
    CloudFunctionObject<CloudType>(pzi),
    writeFile(pzi),
    cellZoneName_(pzi.cellZoneName_),
    cellZoneId_(pzi.cellZoneId_),
    data_(pzi.data_),
    movedParticles_(pzi.movedParticles_),
    maxIDs_(Pstream::nProcs()),
    writerPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleZoneInfo<CloudType>::preEvolve
(
    const typename parcelType::trackingData& td
)
{}


template<class CloudType>
void Foam::ParticleZoneInfo<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    Log_<< this->type() << ":" << nl
        << "    Cell zone                       = " << cellZoneName_ << nl
        << "    Contributions                   = "
        << returnReduce(movedParticles_.size(), sumOp<label>())
        << endl;

    if (!this->writeTime())
    {
        Log_<< endl;
    }

    for (const auto& p : movedParticles_)
    {
        const label id = getParticleID(p);

        if (id == -1)
        {
            // New particle
            data_.append(p);
            maxIDs_[p.origProc] = max(maxIDs_[p.origProc], p.origID);
        }
        else
        {
            // Add to existing particle
            data_[id] += p;
        }
    }

    movedParticles_.clear();

    // Calls write
    CloudFunctionObject<CloudType>::postEvolve(td);
}


template<class CloudType>
bool Foam::ParticleZoneInfo<CloudType>::postMove
(
    parcelType& p,
    const scalar dt,
    const point&,
    const typename parcelType::trackingData& td
)
{
    if (inZone(p.cell()))
    {
        particleInfo newData;
        newData.origID = p.origId();
        newData.origProc = p.origProc();
        newData.position = p.position();
        newData.time0 = this->owner().time().value() + dt;
        newData.age = dt;
        newData.d0 = p.d();
        newData.d = p.d();
        newData.mass0 = p.mass();
        newData.mass = newData.mass0;

        movedParticles_.append(newData);
    }

    return true;
}


template<class CloudType>
void Foam::ParticleZoneInfo<CloudType>::write()
{
    autoPtr<OFstream> osPtr =
        this->newFileAtTime
        (
            "particles",
            this->owner().time().timeOutputValue()
        );

    if (Pstream::parRun())
    {
        // Find number of particles per proc
        labelList allMaxIDs(maxIDs_);
        Pstream::listCombineReduce(allMaxIDs, maxEqOp<label>());

        // Combine into single list
        label n = returnReduce(data_.size(), sumOp<label>());
        DynamicList<particleInfo> globalParticles(n);
        {
            List<List<particleInfo>> procParticles(Pstream::nProcs());
            forAll(procParticles, proci)
            {
                procParticles[proci].resize(allMaxIDs[proci] + 1);
            }

            // Insert into bins for accumulation
            for (const auto& d : data_)
            {
                procParticles[d.origProc][d.origID] = d;
            }

            for (auto& particles : procParticles)
            {
                Pstream::listCombineGather(particles, particleInfoCombineOp());
                for (const auto& p : particles)
                {
                    if (p.origID != -1)
                    {
                        globalParticles.append(p);
                    }
                }
            }
        }

        if (Pstream::master())
        {
            writeWriter(globalParticles);

            auto& os = osPtr();
            writeFileHeader(os);

            label nData = 0;

            for (const auto& p : globalParticles)
            {
                if (p.origID != -1)
                {
                    os << p << endl;
                    ++nData;
                }
            }

            Log_<< "    Number of particles             = " << nData << nl
                << "    Written data to " << os.name() << endl;

            this->setModelProperty("data", globalParticles);
            this->setModelProperty("maxIDs", allMaxIDs);
        }
        else
        {
            // Data only present on master
            this->setModelProperty("data", List<particleInfo>());
            this->setModelProperty("maxIDs", labelList());
        }
    }
    else
    {
        writeWriter(data_);

        auto& os = osPtr();
        writeFileHeader(os);

        for (const auto& p : data_)
        {
            os << p << nl;
        }

        Log_<< "    Number of particles             = " << data_.size() << nl
            << "    Written data to " << os.name() << endl;

        this->setModelProperty("data", data_);
        this->setModelProperty("maxIDs", maxIDs_);
    }

    Log_<< endl;
}


// ************************************************************************* //
