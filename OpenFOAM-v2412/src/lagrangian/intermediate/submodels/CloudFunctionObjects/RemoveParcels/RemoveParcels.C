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

#include "RemoveParcels.H"
#include "fvMesh.H"
#include "faceZone.H"
#include "OFstream.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::RemoveParcels<CloudType>::makeLogFile
(
    const word& zoneName,
    const label zoneI,
    const label nFaces,
    const scalar totArea
)
{
    // Create the output file if not already created
    if (logToFile_)
    {
        DebugInfo<< "Creating output file." << endl;

        if (Pstream::master())
        {
            // Create directory if does not exist
            mkDir(this->writeTimeDir());

            // Open new file at start up
            outputFilePtr_.set
            (
                zoneI,
                new OFstream
                (
                    this->writeTimeDir()/(type() + '_' + zoneName + ".dat")
                )
            );

            outputFilePtr_[zoneI]
                << "# Source    : " << type() << nl
                << "# Face zone : " << zoneName << nl
                << "# Faces     : " << nFaces << nl
                << "# Area      : " << totArea << nl
                << "# Time" << tab << "nParcels" << tab << "mass" << endl;
        }
    }
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::RemoveParcels<CloudType>::postEvolve
(
    const typename parcelType::trackingData& td
)
{
    Log_<< this->modelName() << " output:" << nl;

    const fvMesh& mesh = this->owner().mesh();
    const faceZoneMesh& fzm = mesh.faceZones();

    forAll(faceZoneIDs_, i)
    {
        const word& zoneName = fzm[faceZoneIDs_[i]].name();

        scalar zoneMass = returnReduce(mass_[i], sumOp<scalar>());
        label zoneNParcels = returnReduce(nParcels_[i], sumOp<label>());

        Log_<< "    faceZone " << zoneName
            << ": removed " << zoneNParcels
            << " parcels with mass " << zoneMass
            << nl;
    }

    CloudFunctionObject<CloudType>::postEvolve(td);
}


template<class CloudType>
void Foam::RemoveParcels<CloudType>::write()
{
    const fvMesh& mesh = this->owner().mesh();
    const Time& time = mesh.time();


    List<scalar> allZoneMass(faceZoneIDs_.size(), 0.0);
    List<label> allZoneNParcels(faceZoneIDs_.size(), 0);

    forAll(faceZoneIDs_, i)
    {
        allZoneMass[i] = returnReduce(mass_[i], sumOp<scalar>());
        allZoneNParcels[i] = returnReduce(nParcels_[i], sumOp<label>());

        if (outputFilePtr_.set(i))
        {
            OFstream& os = outputFilePtr_[i];
            os  << time.timeName() << token::TAB
                << allZoneNParcels[i] << token::TAB
                << allZoneMass[i] << endl;
        }
    }

    Log_<< endl;

    if (resetOnWrite_)
    {
        forAll(mass_, i)
        {
            mass_[i] = 0.0;
            nParcels_[i] = 0.0;
        }
    }

    this->setModelProperty("mass", allZoneMass);
    this->setModelProperty("nParcels", allZoneNParcels);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::RemoveParcels<CloudType>::RemoveParcels
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
    faceZoneIDs_(),
    nParcels_(),
    mass_(),
    typeId_(this->coeffDict().template getOrDefault<label>("parcelType", -1)),
    logToFile_(this->coeffDict().getBool("log")),
    resetOnWrite_(this->coeffDict().getBool("resetOnWrite")),
    resetOnStart_(this->coeffDict().getBool("resetOnStart")),
    outputFilePtr_()
{
    const wordList faceZoneNames(this->coeffDict().lookup("faceZones"));

    nParcels_.setSize(faceZoneNames.size(), 0);
    mass_.setSize(faceZoneNames.size(), 0.0);

    if (!resetOnStart_ && Pstream::master())
    {
        this->getModelProperty("mass", mass_);
        this->getModelProperty("nParcels", nParcels_);
    }

    outputFilePtr_.setSize(faceZoneNames.size());

    DynamicList<label> zoneIDs;
    const faceZoneMesh& fzm = owner.mesh().faceZones();
    const surfaceScalarField& magSf = owner.mesh().magSf();
    const polyBoundaryMesh& pbm = owner.mesh().boundaryMesh();

    forAll(faceZoneNames, i)
    {
        const word& zoneName = faceZoneNames[i];
        label zonei = fzm.findZoneID(zoneName);

        if (zonei != -1)
        {
            zoneIDs.append(zonei);
            const faceZone& fz = fzm[zonei];

            label nFaces = returnReduce(fz.size(), sumOp<label>());
            Info<< "        " << zoneName << " faces: " << nFaces << nl;

            scalar totArea = 0.0;
            for (const label facei : fz)
            {
                if (facei < owner.mesh().nInternalFaces())
                {
                    totArea += magSf[facei];
                }
                else
                {
                    const label patchi = pbm.patchID(facei);
                    const polyPatch& pp = pbm[patchi];

                    if
                    (
                        !magSf.boundaryField()[patchi].coupled()
                     || refCast<const coupledPolyPatch>(pp).owner()
                    )
                    {
                        label localFacei = pp.whichFace(facei);
                        totArea += magSf.boundaryField()[patchi][localFacei];
                    }
                }
            }
            totArea = returnReduce(totArea, sumOp<scalar>());

            makeLogFile(zoneName, i, nFaces, totArea);
        }
    }

    faceZoneIDs_.transfer(zoneIDs);
}


template<class CloudType>
Foam::RemoveParcels<CloudType>::RemoveParcels
(
    const RemoveParcels<CloudType>& rpf
)
:
    CloudFunctionObject<CloudType>(rpf),
    faceZoneIDs_(rpf.faceZoneIDs_),
    nParcels_(rpf.nParcels_),
    mass_(rpf.mass_),
    typeId_(rpf.typeId_),
    logToFile_(rpf.logToFile_),
    resetOnWrite_(rpf.resetOnWrite_),
    resetOnStart_(rpf.resetOnStart_),
    outputFilePtr_()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::RemoveParcels<CloudType>::postFace
(
    const parcelType& p,
    const typename parcelType::trackingData& td
)
{
    bool keepParticle = true;

    if ((typeId_ >= 0) && (p.typeId() != typeId_))
    {
        // Not processing this particle type
        return keepParticle;
    }

    if
    (
        this->owner().solution().output()
     || this->owner().solution().transient()
    )
    {
        const faceZoneMesh& fzm = this->owner().mesh().faceZones();

        forAll(faceZoneIDs_, i)
        {
            const faceZone& fz = fzm[faceZoneIDs_[i]];
            if (fz.found(p.face()))
            {
                nParcels_[i]++;
                mass_[i] += p.mass()*p.nParticle();
                keepParticle = false;
                break;
            }
        }
    }

    return keepParticle;
}


// ************************************************************************* //
