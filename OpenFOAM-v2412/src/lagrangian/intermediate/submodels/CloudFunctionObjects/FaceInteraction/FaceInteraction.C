/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021-2022 OpenCFD Ltd.
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

#include "FaceInteraction.H"
#include "faceZoneMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<class CloudType>
const Foam::Enum<typename Foam::FaceInteraction<CloudType>::interactionType>
Foam::FaceInteraction<CloudType>::interactionTypeNames_
({
    { interactionType::STICK, "stick" },
    { interactionType::ESCAPE, "escape" },
    { interactionType::REBOUND, "rebound" },
});


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
bool Foam::FaceInteraction<CloudType>::processParticle
(
    const parcelType& p,
    const label localZonei
)
{
    if (!faceZoneBBs_[localZonei].contains(p.position()))
    {
        // Quick reject if the particle is not in the face zone bound box
        return false;
    }

    if ((p.d() > dMax_) || (p.d() < dMin_))
    {
        return false;
    }

    return true;
}


template<class CloudType>
void Foam::FaceInteraction<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const faceZoneMesh& fzm = mesh.faceZones();

    Log_<< type() << " output:" << nl;

    // Retrieve any stored data
    const label nZones = faceZoneIDs_.size();
    labelList npe0(nZones, Zero);
    labelList nps0(nZones, Zero);
    labelList npr0(nZones, Zero);

    this->getModelProperty("nEscape", npe0);
    this->getModelProperty("nStick", nps0);
    this->getModelProperty("nRebound", npr0);

    // Accumulate current data
    labelList npe(returnReduce(nEscapeParticles_, sumOp<labelList>()));
    labelList nps(returnReduce(nStickParticles_, sumOp<labelList>()));
    labelList npr(returnReduce(nReboundParticles_, sumOp<labelList>()));
    forAll(npe, i)
    {
        npe[i] = npe[i] + npe0[i];
        nps[i] = nps[i] + nps0[i];
        npr[i] = npr[i] + npr0[i];
    }


    // Write to log/file
    forAll(faceZoneIDs_, i)
    {
        const label zonei = faceZoneIDs_[i];
        Log_<< "    Zone : " << fzm[zonei].name() << nl
            << "        Escape  : " << npe[i] << nl
            << "        Stick   : " << nps[i] << nl
            << "        Rebound : " << npr[i] << nl;

        if (this->writeToFile())
        {
            auto& os = filePtrs_[i];

            writeCurrentTime(os);

            // Note - writing as scalar for better formatting
            os  << tab << scalar(npe[i])
                << tab << scalar(nps[i])
                << tab << scalar(npr[i])
                << endl;
        }
    }
    Log_<< endl;

    // Set restart data
    this->setModelProperty("nEscape", npe);
    this->setModelProperty("nStick", nps);
    this->setModelProperty("nRebound", npr);

    nEscapeParticles_ = Zero;
    nStickParticles_ = Zero;
    nReboundParticles_ = Zero;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::FaceInteraction<CloudType>::FaceInteraction
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
        typeName,
        this->coeffDict()
    ),
    faceZoneIDs_(),
    faceZoneBBs_(),
    faceZoneInteraction_(),
    filePtrs_(),
    nEscapeParticles_(),
    nStickParticles_(),
    nReboundParticles_(),
    dMin_(this->coeffDict().getOrDefault("dMin", -GREAT)),
    dMax_(this->coeffDict().getOrDefault("dMax", GREAT))
{
    const List<Tuple2<word, word>> nameAndInteraction
    (
        this->coeffDict().lookup("faceZones")
    );

    filePtrs_.setSize(nameAndInteraction.size());

    faceZoneBBs_.setSize(nameAndInteraction.size());
    faceZoneInteraction_.setSize(nameAndInteraction.size());

    DynamicList<label> zoneIDs;
    const faceZoneMesh& fzm = owner.mesh().faceZones();
    const auto& faces = this->owner().mesh().faces();
    const auto& points = this->owner().mesh().points();

    label nZone = 0;
    for (const auto& zoneInfo : nameAndInteraction)
    {
        const word& zoneName = zoneInfo.first();
        const label zonei = fzm.findZoneID(zoneName);

        if (zonei != -1)
        {
            zoneIDs.append(zonei);
            const faceZone& fz = fzm[zonei];

            const label nFaces = returnReduce(fz.size(), sumOp<label>());
            Info<< "        " << zoneName << " faces: " << nFaces << nl;

            // Cache the particle action for this faceZone
            faceZoneInteraction_[nZone] =
                interactionTypeNames_[zoneInfo.second()];

            // Cache faceZone bound box
            auto& bb = faceZoneBBs_[nZone];
            for (const label facei : fz)
            {
                for (const label fpi : faces[facei])
                {
                    bb.add(points[fpi]);
                }
            }
            reduce(bb.min(), minOp<point>());
            reduce(bb.max(), maxOp<point>());
            bb.inflate(0.05);

            // Create output file for zone
            if (this->writeToFile())
            {
                filePtrs_.set
                (
                    nZone,
                    this->newFileAtStartTime(modelName + '_' + zoneName)
                );

                writeHeaderValue(filePtrs_[nZone], "Source", type());
                writeHeaderValue(filePtrs_[nZone], "Face zone", zoneName);
                writeHeaderValue(filePtrs_[nZone], "Faces", nFaces);
                writeCommented(filePtrs_[nZone], "Time");
                writeTabbed(filePtrs_[nZone], "Escape");
                writeTabbed(filePtrs_[nZone], "Stick");
                writeTabbed(filePtrs_[nZone], "Rebound");
                filePtrs_[nZone] << nl;
            }
            ++nZone;
        }
        else
        {
            WarningInFunction
                << "Unable to find faceZone " << zoneName
                << " - removing" << endl;
        }
    }
    faceZoneIDs_.transfer(zoneIDs);

    filePtrs_.setSize(nZone);
    faceZoneBBs_.setSize(nZone);
    faceZoneInteraction_.setSize(nZone);
    nEscapeParticles_.setSize(nZone, Zero);
    nStickParticles_.setSize(nZone, Zero);
    nReboundParticles_.setSize(nZone, Zero);
}


template<class CloudType>
Foam::FaceInteraction<CloudType>::FaceInteraction
(
    const FaceInteraction<CloudType>& pfi
)
:
    CloudFunctionObject<CloudType>(pfi),
    functionObjects::writeFile(pfi),
    faceZoneIDs_(pfi.faceZoneIDs_),
    faceZoneBBs_(pfi.faceZoneBBs_),
    filePtrs_(),
    nEscapeParticles_(pfi.nEscapeParticles_),
    nStickParticles_(pfi.nStickParticles_),
    nReboundParticles_(pfi.nReboundParticles_),
    dMin_(pfi.dMin_),
    dMax_(pfi.dMax_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::FaceInteraction<CloudType>::postFace
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    bool keepParticle = true;

    const auto& fzm = this->owner().mesh().faceZones();

    forAll(faceZoneIDs_, i)
    {
        if (!processParticle(p, i))
        {
            continue;
        }

        const label zonei = faceZoneIDs_[i];

        const label localFacei = fzm[zonei].find(p.face());

        if (localFacei != -1)
        {
            const label facei = fzm[zonei][localFacei];

            switch (faceZoneInteraction_[i])
            {
                case interactionType::ESCAPE:
                {
                    keepParticle = false;
                    ++nEscapeParticles_[i];
                    break;
                }
                case interactionType::STICK:
                {
                    auto& pRef = const_cast<parcelType&>(p);
                    pRef.active(false);
                    pRef.U() = Zero;
                    ++nStickParticles_[i];
                    break;
                }
                case interactionType::REBOUND:
                {
                    const face& f = this->owner().mesh().faces()[facei];
                    const auto n = f.unitNormal(this->owner().mesh().points());

                    auto& pRef = const_cast<parcelType&>(p);
                    pRef.U() -= 2*n*(pRef.U() & n);
                    ++nReboundParticles_[i];
                    break;
                }
                default:
                {
                    FatalErrorInFunction
                        << "Unhandled enumeration "
                        << interactionTypeNames_[faceZoneInteraction_[i]]
                        << abort(FatalError);
                }
            }
        }
    }

    return keepParticle;
}


// ************************************************************************* //
